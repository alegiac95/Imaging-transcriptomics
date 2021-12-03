import numpy as np
from pathlib import Path
from scipy.stats import zscore, spearmanr
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from netneurotools import freesurfer, stats
from pyls import pls_regression
import pandas as pd
from .inputs import load_gene_expression, \
    extract_average, \
    read_scan, \
    load_gene_labels
from statsmodels.stats.multitest import multipletests
from scipy.stats import norm
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from functools import partial
from itertools import product
import gseapy


np.random.seed(1234)


class ImagingTranscriptomics:
    # --------- INITIALIZATION --------- #
    def __init__(self,
                 scan_data,
                 regions="cort+sub",
                 method="pls",
                 **kwargs):
        """ImagingTranscriptomics class for imaging transcriptomics analysis.

        :param np.array scan_data: imaging scan data.
        :param str regions: regions to be used for analysis.
        :param str method: method to run the analysis, can be either "pls"
        for pls regression or "corr" cor simple correlation analysis.
        :param kwargs: additional arguments for the method. This include:
            * "n_components": number of components for pls regression.
            * "var": variance explained for pls regression.
            * "n_permutations": number of permutations for permutation test.

        """
        if regions == "cort+sub" or regions == "all":
            assert scan_data.shape == (41,)
            self.zscore_data = zscore(scan_data, axis=0, ddof=1)
            self._cortical = self.zscore_data[:34]
            self._subcortical = self.zscore_data[34:]
        elif regions == "cort":
            assert scan_data.shape == (34,) or scan_data.shape == (41,)
            self.zscore_data = zscore(scan_data, axis=0, ddof=1)
            self._cortical = self.zscore_data if scan_data.shape == (34,) else\
                self.zscore_data[34:]
            self._subcortical = None
        self._regions = regions
        self.scan_data = scan_data
        if method not in ["pls", "corr"]:
            raise ValueError(
                "The method must be either pls or corr."
                "Please choose either pls or corr to run the analysis."
            )
        else:
            self._method = method
        if self._method == "pls":
            if "n_components" not in kwargs and "var" not in kwargs:
                raise ValueError("You must specify either the variance or "
                                 "the number of components for pls regression")
            else:
                self._analysis = PLSAnalysis(kwargs.get("n_components"),
                                             kwargs.get("var"))
        elif self._method == "corr":
            self._analysis = CorrAnalysis()
        self.gene_expression = load_gene_expression(self._regions)
        self.gene_labels = load_gene_labels()
        self._permutations = None
        self._permutation_ind = None

    @classmethod
    def from_scan(cls, scan_path, method="pls", regions="cort+sub", **kwargs):
        """Initialise an ImagingTranscriptomics object from a NIfTI scan.
        The extracted data corresponds to the average of the ROIs in the DK
        atlas.
        """
        scan_data = extract_average(read_scan(scan_path))
        return cls(scan_data, method=method, regions=regions, **kwargs)

    @classmethod
    def from_file(cls, file_path, method="pls", regions="cort+sub", **kwargs):
        """Initialise an ImagingTranscriptomics object from a text file.
        The file should contain a column with the data you want to use for
        the analysis.
        """
        scan_data = np.loadtxt(file_path)
        return cls(scan_data, method=method, regions=regions, **kwargs)

    # --------- PROPERTIES --------- #
    @property
    def method(self):
        return self._method

    @property
    def gene_results(self):
        return self._analysis.gene_results

    # --------- METHODS --------- #
    def permute_data(self, n_permutations=1000):
        """Permute the imaging data maintaining spatial autocorrelation for
        the cortical regions. The permutation is done using the 
        netneurotools Python package.

        :param int n_permutations: number of permutations.
        """
        _permuted = np.zeros((self.zscore_data.shape[0], n_permutations))
        _perm_indexes = np.zeros((self.zscore_data.shape[0], n_permutations))
        # Calculate the permutations on the subcortical regions.
        if self._subcortical is not None:
            sub_permuted = np.zeros((self._subcortical.shape[0],
                                     n_permutations))
            for i in range(n_permutations):
                sub_resample = np.random.choice(7, size=7)
                _perm_indexes[34:, i] = sub_resample + 34  # keep into
                # account the shift of the subcortical given by the cortical
                # regions.
                sub_permuted[:, i] = self._subcortical[sub_resample]
            _permuted[34:, :] = sub_permuted
        # Cortical
        # Annotation file for the Desikan-Killiany atlas in fs5
        annot_lh = Path(__file__).resolve().parent / "data/fsa5_lh_aparc.annot"
        annot_rh = Path(__file__).resolve().parent / "data/fsa5_rh_aparc.annot"
        # Get the parcel centroids of the Desikan-Killiany atlas
        parcel_centroids, parcel_hemi = freesurfer.find_parcel_centroids(
            lhannot=annot_lh,
            rhannot=annot_rh,
            version="fsaverage5",
            surf="sphere",
            method="surface",
        )
        # Mask the results to have only the left hemisphere
        left_hemi_mask = parcel_hemi == 0
        parcel_centroids, parcel_hemi = (
            parcel_centroids[left_hemi_mask],
            parcel_hemi[left_hemi_mask],
        )
        # Get the spin samples
        spins = stats.gen_spinsamples(
            parcel_centroids, parcel_hemi,
            n_rotate=n_permutations,
            method="vasa",
            seed=1234
        )
        cort_permuted = np.array(self._cortical[spins]).reshape(34,
                                                                n_permutations)
        _perm_indexes[:34, :] = spins
        _permuted[0:34, :] = cort_permuted
        self._permutations = _permuted
        self._permutation_ind = _perm_indexes

    # --------- RUN ANALYSIS --------- #
    def run(self, outdir=None, gsea=True, gene_set="lake", noplot=False):
        """Method to run the imaging transcriptomics analysis.

        :param str outdir: path to the output directory, if not provided the
        results will be saved in the current directory.
        :param bool gsea: if True, run the GSEA analysis, if False the GSEA
        analysis is skipped.
        :param str gene_set: gene set to use for the GSEA analysis.
        :param bool noplot: if True, do not plot the results.
        """
        self.permute_data()
        if self._method == "corr":
            if self._regions == "cort":
                _d = self._cortical
                if self._permutations.shape[0] == 41:
                    _d_perm = self._permutations[0:34, :]
                else:
                    _d_perm = self._permutations
            elif self._regions == "cort+sub" or self._regions == "all":
                _d = self.zscore_data
                _d_perm = self._permutations
            self._analysis.bootstrap_correlation(_d, _d_perm,
                                                 self.gene_expression)
            if gsea:
                self._analysis.gsea(gene_set=gene_set, outdir=outdir)

        pass  # TODO: make final run() method.


# --------- GENE ANALYSIS --------- #
class GeneResults:
    def __init__(self, method, **kwargs):
        """Initialize the results of the analysis. Depending on the method
        used, the results will have underlying result classes, which account
        for the different analysis methods.

        :param str method: the method used for the analysis.
        :param kwargs: Additional parameters, for the initialisation. If the
        method is "pls" among the kwargs you *MUST* specify the number of
        components used, for the initialisation of the pls class.
        """
        self.method = method
        if self.method == "pls":
            self.results = PLSGenes(kwargs.get("n_components"))
        elif self.method == "corr":
            self.results = CorrGenes()
        else:
            raise ValueError(f"The method {method} is not supported.")

    @property
    def n_genes(self):
        return self.results.n_genes

    @property
    def genes(self):
        if isinstance(self.results, PLSGenes):
            pass
        elif isinstance(self.results, CorrGenes):
            return self.results.genes

    @property
    def scores(self):
        if isinstance(self.results, PLSGenes):
            pass
        elif isinstance(self.results, CorrGenes):
            return self.results.corr

    @property
    def boot(self):
        if isinstance(self.results, PLSGenes):
            pass
        elif isinstance(self.results, CorrGenes):
            return self.results.boot_corr

    @property
    def pvals(self):
        if isinstance(self.results, PLSGenes):
            pass
        elif isinstance(self.results, CorrGenes):
            return self.results.pval

    @property
    def pvals_corr(self):
        if isinstance(self.results, PLSGenes):
            pass
        elif isinstance(self.results, CorrGenes):
            return self.results.pval_corr


# --------- PLS ANALYSIS --------- #
class PLSAnalysis:
    """Class for performing PLS regression on the imaging data.
    This class contains fields for
    """
    def __init__(self, n_components: int, var: float):
        self.n_components = n_components
        self.var = self.check_var(var)
        if self.var is None and self.n_components is None:
            raise ValueError("You must specify either the variance or "
                             "the number of components for pls regression")
        self.components_var = np.zeros(15)
        self._p_val = np.zeros(self.n_components)
        self._r2 = np.zeros(self.n_components)
        self.gene_results = GeneResults("pls", n_components=self.n_components)

    @staticmethod
    def check_var(var: float):
        """ Check if the variance is between 0 and 1.

        :param float var: variance.
        """
        if var < 0:
            raise ValueError("The variance must be a positive number.")
        elif var > 1:
            raise ValueError("The variance must be a number between 0 and 1.")
        return var

    def set_coef(self, data, gene_exp):  # TODO: need to check
        """Set the coefficients for the PLS regression. The function will
        estimate the variance or the number of components depending on the
        non missing parameter through a PLS regression with 15 components.

        :param data: imaging data.
        :param gene_exp: gene expression data.
        """
        res = pls_regression(data, gene_exp,
                             n_components=15,
                             n_perm=0,
                             n_boot=0)
        explained_var = res.get("varexp")
        if self.var is None:
            self.var = np.cumsum(explained_var)[self.n_components-1]
        if self.n_components is None:
            dim = 1
            cumulative_var = np.cumsum(explained_var)
            while cumulative_var[dim - 1] < self.var:
                dim += 1
            self.n_components = dim
        return explained_var

    @property
    def p_val(self):
        """P-value of the pls regression given by the bootstrapping on the
        PLS components."""
        return self._p_val

    @property
    def r2(self):
        """R2 of the pls regression given by the bootstrapping on the PLS
        components"""
        return self._r2

    @p_val.setter
    def p_val(self, p_val):
        self._p_val = p_val

    @r2.setter
    def r2(self, r2):
        self._r2 = r2

    def boot_pls(self, imaging_data, permuted_imaging, gene_exp):
        """Bootstrapping on the PLS components.

        :param imaging_data: imaging data. Allows the user to specify the
        data to use (e.g., with only cortical regions this can be only the
        cortical vector, other wise the whole data).
        :param permuted_imaging: imaging data permuted. Allows the user to
        specify the data to use (e.g., with only cortical regions this can be
        only the cortical vector, other wise the whole data).
        :param gene_exp: gene expression data.
        """
        # Iterate over the components
        for component in range(1, self.n_components + 1):
            _res = pls_regression(imaging_data, gene_exp,
                                  n_components=component, n_perm=0, n_boot=0)
            _exp_var = _res.get("varexp")
            _temp = _exp_var.cumsum(axis=0)
            _R = _temp[component - 1]
            _R_sq = np.zeros(1000)
            # Iterate over the permutations
            for i in tqdm(range(1000),
                          desc=f"Bootstrapping on PLS component {component}",
                          unit=" iterations"):
                y_data = permuted_imaging[:, i].reshape(
                    imaging_data.shape[0], 1)
                _res = pls_regression(y_data, gene_exp, n_components=component,
                                      n_perm=0, n_boot=0)
                _exp_var = _res.get("varexp")
                _R_sq[i] = _exp_var.cumsum(axis=0)[component - 1]
            # Set the results
            self.r2[component - 1] = _R
            self.p_val[component - 1] = np.sum(_R_sq >= _R) / 1000


# --------- PLS GENES --------- #
class PLSGenes:
    def __init__(self, n_components):
        """ Initialize the results of the PLS analysis. The result will
        include both the permuted and the original results. The class
        contains two subclasses, one for the original results and one for
        the bootrapped results.

        :param int n_components: number of components used for the analysis.
        """
        self.n_genes = 15633
        self.n_components = n_components
        self.orig = OrigPLS(n_components, self.n_genes)
        self.boot = BootPLS(n_components, self.n_genes)


# --------- ORIG PLS --------- #
class OrigPLS:
    def __init__(self, n_components, n_genes):
        """ Initialize the original results of the PLS analysis. The class
        contains the fields corresponding to the number of components used,
        the weights of the pls for each gene ordered in descending order,
        the index where the original genes and the zscore of the weights.

        :param int n_components: number of components used.
        :param int n_genes: number of genes.
        """
        self.n_components = n_components
        self.weights = np.zeros((n_components, n_genes))
        self.genes = np.zeros((n_components, n_genes))
        self.index = np.zeros((n_components, n_genes))
        self.zscored = np.zeros((n_components, n_genes))


# --------- BOOT PLS --------- #
class BootPLS:
    def __init__(self, n_components, n_genes, n_iter=1000):
        """Initialise a class to store the results of the bootstrapping of
        the genes.

        All the initialised fields are stored as numpy arrays with the
        number of rows corresponding to the number of components and the
        columns corresponding to the number of genes. The weights field has
        an additional 3rd dimension corresponding to the number of
        iterations (the default number is 1000).
        The fields are:

        * weights (n_components, n_genes, n_iter): the weights of the genes
        for each component, for each iteration.

        * genes (n_components, n_genes, n_iter): the genes that correspond
        to the most contributing genes for each component.

        * index (n_components, n_genes, n_iter): the index of the genes
        compared to the original list of gene labels.

        * std: the standard deviation of the weights for each component,
        calculated from the bootstrapped weights.

        * zscored (n_components, n_genes, n_iter): the z-scored weights.

        * pval (n_components, n_genes): the p-value of the z-scored gene
        wights.

        * pval_corr (n_components, n_genes): the p-value of the correlation
        corrected for multiple comparisons using the Benjamini-Hochberg method.

        :param int n_components: number of components used for the analysis.
        :param int n_genes: number of genes used for the analysis.
        :param int n_iter: number of iterations used for the bootstrapping,
        the default is 1000.
        """
        self.n_components = n_components
        self.weights = np.zeros((n_components, n_genes, n_iter))
        self.genes = np.zeros((n_components, n_genes))
        self.std = np.zeros((n_components, n_genes))
        self._z_score = np.zeros((n_components, n_genes))
        self.pval = np.zeros((n_components, n_genes))
        self.pval_corr = np.zeros((n_components, n_genes))

    @property
    def z_score(self):
        """The z-scored weights.
        """
        return self._z_score

    def compute(self, weights, genes):
        """ Compute the p-values of the z-scored weights.
        The compute function calculates the zscores based on the original
        weights, orders the weights in descending order and then calculates
        the p-values and corrects for multiple comparisons using the
        Benjamini-Hochberg method.

        :param numpy.ndarray weights: the weights of the genes for each of
        the components. The array should have shape (n_components, n_genes).
        :param numpy.ndarray genes: the labels in the original order.
        """
        for component in range(self.n_components):
            # Calculate the standard deviation of all the wights
            self.std[component, :] = self.weights[component, :, :].std(
                axis=0, ddof=1)
            # Calculate the zscore and store it sorted in descending order
            _z = weights[component, :] / self.std[component, :]
            self.z_score[component, :] = np.sort(
                _z, axis=0, kind='mergesort')[::-1]
            _index = np.argsort(_z, axis=0, kind='mergesort')[::-1]
            # Reorder the genes according to the zscore
            self.genes[component, :] = genes[_index]
            # Calculage pvalue and pvalue corrected
            _p_val = norm.sf(abs(self.z_score[component, :]))
            self.pval[component, :] = _p_val
            _, _p_corr, _, _ = multipletests(_p_val, method='fdr_bh')
            self.pval_corr[component, :] = _p_corr


# --------- CORR ANALYSIS --------- #
def _spearman_op(idx, permuted, y):
    """Wrapper for the spearman correlation function, to be used for
    parallel computation."""
    return spearmanr(permuted[:, idx[0]], y[:, idx[1]])[0]


class CorrAnalysis:
    """Class to store all the results of the correlation analysis.
    The only field contained is the a GeneResults object, storing all
    information about the correlation analysis.
    """
    def __init__(self):
        self.gene_results = GeneResults("corr")
        self.gsea_res = None

    # --------- COMPUTE FUNCTIONS --------- #
    def bootstrap_correlation(self, imaging_data, permuted_imaging, gene_exp):
        """Perform bootstrapping on the correlation.

        The function first calculates the correlation between the imaging
        vector and each of the genes. Then, it performs 1000 bootstrapping
        iterations of the same correaltion only using the permuted imaging
        data. The function then calls the methods of the GeneResults class
        to order the genes and calculate the p-values.

        :param np.ndarray imaging_data: the imaging data. This vector
        represents the original data to be correlated, and can have either
        length 41 or 34, depending on whether the subcortical regions want
        to be included in the analysis.
        :param np.ndarray permuted_imaging: the imaging data permuted. This
        is a matrix with shape (n_imaging, 1000), where `n_imaging`
        represents the length of the imaging vector.
        :param np.ndarray gene_exp: the gene expression data.
        """
        assert isinstance(self.gene_results.results, CorrGenes)
        for i in range(self.gene_results.n_genes):
            self.gene_results.results.corr[0, i], _ = spearmanr(
                imaging_data, gene_exp[:, i])
        pool = Pool(cpu_count()-1)
        _ind = product(range(1000), range(self.gene_results.n_genes))

        print("Performing bootstrapping...")
        for ind, res in enumerate(pool.imap(partial(_spearman_op,
                                                    permuted=permuted_imaging,
                                                    y=gene_exp),
                                            _ind, chunksize=1000)):
            self.gene_results.results.boot_corr.flat[ind] = res
        self.gene_results.results.sort_genes()
        self.gene_results.results.compute_pval()

    def gsea(self, gene_set="lake", outdir=None):
        """Perform GSEA on the correlation."""
        assert isinstance(self.gene_results.results, CorrGenes)
        gene_set = get_geneset(gene_set)
        # prepare the gene_list as a list of strings
        gene_list = [
            str(gene) for gene in self.gene_results.results.genes[1, :]
        ]
        # perform the GSEA on orignal results
        rnk = pd.DataFrame(list(
                zip(gene_list, self.gene_results.results.corr[1, :])))
        gsea_results = gseapy.prerank(rnk, gene_set,
                                      outdir=None,
                                      permutation_num=100,
                                      seed=1234)
        _origin_es = gsea_results.res2d.es.to_numpy()
        _boot_es = np.zeros((_origin_es.shape[0], 1000))
        # perform the GSEA on the permutations
        for i in range(1000):
            rnk = pd.DataFrame(list(
                zip(gene_list, self.gene_results.results.boot_corr[:, i])
            ))
            _gsea_res = gseapy.prerank(rnk, gene_set,
                                       permutation_num=1,
                                       no_plot=True,
                                       outdir=None,
                                       seed=1234)
            _boot_es[:, i] = _gsea_res.res2d.es.to_numpy()
        # calculate the p-value
        _p_val = np.zeros((_origin_es.shape[0],))
        for i in range(_origin_es.shape[0]):
            _p_val[i] = np.sum(_boot_es[i, :] >= _origin_es[i]) / 1000
        # calculate the p-value corrected
        _, _p_corr, _, _ = multipletests(_p_val, method='fdr_bh')
        self.gsea_res = gsea_results  # TODO: use the _p_Val and _pval_corr
        # and substitute in the results, so that the results are accurate.
        # TODO: make the output and clean the .csv file
        if outdir is not None:
            outdir = Path(outdir)
            assert outdir.exists()

    # --------- SAVE FUNCTIONS --------- #
    def save_results(self, outdir):
        """Save the results to a file."""
        outdir = Path(outdir)
        assert outdir.exists()
        # Create the data to save
        data_to_save = zip(self.gene_results.results.genes,
                           self.gene_results.results.corr,
                           self.gene_results.results.pval,
                           self.gene_results.results.pval_corr)
        # Save the data
        pd.DataFrame(
            data_to_save, columns=["Gene", "Corr", "Pval", "Pval_corr"]
        ).to_csv(outdir / "corr_genes.tsv", index=False, sep='\t')


def get_geneset(gene_set: str):  # TODO: make a real thing
    """Returns the path to the geneset file, if it is not one of those
    defined in gseapy.

    :param str gene_set: the name of the geneset. If the geneset is "lake"
    or "pooled", the genesets are provided with the package data, otherwise
    the gene sets are from the gseapy package.
    """
    if gene_set.lower() == "lake":
        return str(Path(__file__).parent.parent / "data" / "geneset_LAKE.gmt")
    elif gene_set.lower() == "pooled":
        return str(Path(__file__).parent.parent / "data" /
                   "geneset_Pooled.gmt")
    else:
        return gene_set


# --------- CORRELATION  GENES  --------- #
class CorrGenes:
    """Class that stores the gene results of the correlation analysis. It
    has the following fields:

    * boot_corr: the bootstrapped results of the correlation analysis.
    * corr: the original results of the correlation analysis.
    * boot_corr: the bootstrapped results of the correlation analysis.
    * genes: the gene list used for the analysis.
    * pval: the p-value of the correlation.
    * pval_corr: the p-value of the correlation corrected for multiple
    comparisons using the Benjamini-Hochberg method.
    """
    def __init__(self, n_iter=1000):
        """Initialise the class.

        :param int n_iter: number of iterations used for the bootstrapping,
        """
        self.n_genes = 15633
        self._n_iter = n_iter
        self.boot_corr = np.zeros((self.n_genes, self._n_iter))
        self.corr = np.zeros((1, self.n_genes))
        self.genes = np.zeros((1, self.n_genes))
        self.pval = np.zeros((1, self.n_genes))
        self.pval_corr = np.zeros((1, self.n_genes))
        self._index = None

    def compute_pval(self):
        """Compute the p-values, and its fdr correction, of the correlation,
        from the list of
        bootstrapped correlations.
        """
        # This calculation assumes that the order of the genes is the same
        # in both the original and the bootstrapped list. IF one is ordered,
        # make sure the order of the other is the same.
        for i in range(self.n_genes):
            self.pval[0, i] = np.sum(
                self.boot_corr[i, :] >= self.corr[0, i]) / self._n_iter
        _, p_corr, _, _ = multipletests(self.pval, method='fdr_bh',
                                        is_sorted=False)
        self.pval_corr[0, :] = p_corr

    @property
    def is_sorted(self):
        """Check if the list of genes is sorted."""
        return False if self._index is None else True

    def sort_genes(self):
        """Order the genes in the list of genes. Both the order of the
        order of the bootstrapped genes are ordered.
        """
        self._index = np.argsort(self.corr, axis=1, kind='mergesort')[::-1]
        self.corr[0, :] = self.corr[0, self._index]
        self.genes[0, :] = self.genes[0, self._index]
        self.pval[0, :] = self.pval[0, self._index]
        self.pval_corr[0, :] = self.pval_corr[0, self._index]
        for i in range(self._n_iter):
            self.boot_corr[:, i] = self.boot_corr[:, i][self._index]
