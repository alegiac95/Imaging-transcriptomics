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
from collections import OrderedDict
import gseapy
from gseapy.plot import gseaplot


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
        :param str regions: regions to be used for analysis. These can be
        "cort+sub" (or "all") which will perform the analysis on the
        cortical and subcortical regions, or "cort" which will only perform
        the analysis on the cortical regions.
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
                self.zscore_data[:34]
            self._subcortical = None
        self._regions = regions
        self.scan_data = scan_data
        self.gene_expression = load_gene_expression(self._regions)
        self.gene_labels = load_gene_labels()
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
                if self._regions == "all" or self._regions == "cort+sub":
                    self._analysis = PLSAnalysis(self.zscore_data,
                                                 self.gene_expression,
                                                 kwargs.get("n_components"),
                                                 kwargs.get("var"))
                else:
                    self._analysis = PLSAnalysis(self._cortical,
                                                 self.gene_expression,
                                                 kwargs.get("n_components"),
                                                 kwargs.get("var"))
        elif self._method == "corr":
            self._analysis = CorrAnalysis()
        self._permutations = None
        self._permutation_ind = None

    @classmethod
    def from_scan(cls, scan_path, method="pls", regions="cort+sub", **kwargs):
        """Initialise an ImagingTranscriptomics object from a NIfTI scan.
        The extracted data corresponds to the average of the ROIs in the DK
        atlas.

        :param str scan_path: path to the NIfTI scan.
        :param str method: method to run the analysis, can be either "pls"
        for pls regression or "corr" cor simple correlation analysis.
        :param str regions: regions to be used for analysis. These can be
        "cort+sub" (or "all") which will perform the analysis on the
        cortical and subcortical regions, or "cort" which will only perform
        the analysis on the cortical regions.
        :param kwargs: additional arguments for the method. This include:
            * "n_components": number of components for pls regression.
            * "var": variance explained for pls regression.
            * "n_permutations": number of permutations for permutation test.
        :return: ImagingTranscriptomics object.
        """
        scan_data = extract_average(read_scan(scan_path))
        return cls(scan_data, method=method, regions=regions, **kwargs)

    @classmethod
    def from_file(cls, file_path, method="pls", regions="cort+sub", **kwargs):
        """Initialise an ImagingTranscriptomics object from a text file.
        The file should contain a column with the data you want to use for
        the analysis.

        :param str file_path: path to the text file. The text file should
        contain a column with the data you want to use for the analysis.
        :param str method: method to run the analysis, can be either "pls"
        for pls regression or "corr" cor simple correlation analysis.
        :param str regions: regions to be used for analysis. These can be
        "cort+sub" (or "all") which will perform the analysis on the
        cortical and subcortical regions, or "cort" which will only perform
        the analysis on the cortical regions.
        :param kwargs: additional arguments for the method. This include:
            * "n_components": number of components for pls regression.
            * "var": variance explained for pls regression.
            * "n_permutations": number of permutations for permutation test.
        :return: ImagingTranscriptomics object.
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
        _perm_indexes = np.zeros((self.zscore_data.shape[0],
                                  n_permutations), dtype=np.int32)
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

    @staticmethod
    def _make_output_dir(output_dir):
        outdir = Path(output_dir) / "Imt_"
        outdir.mkdir(exist_ok=True)
        return outdir

    # --------- RUN ANALYSIS --------- #
    def run(self, outdir=None, gsea=True, gene_set="lake"):
        """Method to run the imaging transcriptomics analysis.

        :param str outdir: path to the output directory, if not provided the
        results will be saved in the current directory.
        :param bool gsea: if True, run the GSEA analysis, if False the GSEA
        analysis is skipped.
        :param str gene_set: gene set to use for the GSEA analysis.
        """
        self.permute_data()
        if outdir is not None:
            outdir = self._make_output_dir(outdir)
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
                                                 self.gene_expression,
                                                 self.gene_labels)
            if gsea:
                self._analysis.gsea(gene_set=gene_set, outdir=outdir)
            self._analysis.save_results(outdir=outdir)
        elif self._method == "pls":
            if self._regions == "cort":
                _d = self._cortical
                if self._permutations.shape[0] == 41:
                    _d_perm = self._permutations[0:34, :]
                else:
                    _d_perm = self._permutations
            elif self._regions == "cort+sub" or self._regions == "all":
                _d = self.zscore_data
                _d_perm = self._permutations
            assert isinstance(self._analysis, PLSAnalysis)
            self._ananlysis.boot_pls(_d, _d_perm, self.gene_expression)


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
            return self.results.orig.genes
        elif isinstance(self.results, CorrGenes):
            return self.results.genes

    @property
    def scores(self):
        if isinstance(self.results, PLSGenes):
            return self.results.orig.weights
        elif isinstance(self.results, CorrGenes):
            return self.results.corr

    @property
    def boot(self):
        if isinstance(self.results, PLSGenes):
            return self.results.boot.weights
        elif isinstance(self.results, CorrGenes):
            return self.results.boot_corr

    @property
    def pvals(self):
        if isinstance(self.results, PLSGenes):
            return self.results.boot.pval
        elif isinstance(self.results, CorrGenes):
            return self.results.pval

    @property
    def pvals_corr(self):
        if isinstance(self.results, PLSGenes):
            return self.results.boot.pval_corr
        elif isinstance(self.results, CorrGenes):
            return self.results.pval_corr


# --------- PLS ANALYSIS --------- #
class PLSAnalysis:
    """Class for performing PLS regression on the imaging data.
    This class contains fields for
    """
    def __init__(self, imaging_data, gene_exp,  n_components: int, var: float):
        self.var, self.n_components, self.components_var = self.set_coef(
            imaging_data, gene_exp, n_components=n_components, var=var)
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

    @staticmethod
    def set_coef(data, gene_exp, var=None, n_components=None):
        """Set the coefficients for the PLS regression. The function will
        estimate the variance or the number of components depending on the
        non missing parameter through a PLS regression with 15 components.

        :param data: imaging data.
        :param gene_exp: gene expression data.
        :param var: variance.
        :param n_components: number of components.
        """
        res = pls_regression(gene_exp, data.reshape(data.shape[0], 1),
                             n_components=15,
                             n_perm=0,
                             n_boot=0)
        explained_var = res.get("varexp")
        if var is None:
            var = np.cumsum(explained_var)[n_components-1]
        if n_components is None:
            dim = 1
            cumulative_var = np.cumsum(explained_var)
            while cumulative_var[dim - 1] < var:
                dim += 1
            n_components = dim
        return var, n_components, explained_var

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
            _res = pls_regression(gene_exp, imaging_data.reshape(
                imaging_data.shape[0], 1),
                                  n_components=component,
                                  n_perm=0, n_boot=0)
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
                _res = pls_regression(gene_exp, y_data, n_components=component,
                                      n_perm=0, n_boot=0)
                _exp_var = _res.get("varexp")
                _R_sq[i] = _exp_var.cumsum(axis=0)[component - 1]
            # Set the results
            self.r2[component - 1] = _R
            self.p_val[component - 1] = np.sum(_R_sq >= _R) / 1000
        return


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

    def boot_genes(self, imaging_data, permuted_imaging, perm_indexes,
                   scan_data, gene_exp, gene_labels):
        """Bootstrapping on the PLS components.

        :param imaging_data: imaging data. Allows the user to specify the
        data to use (e.g., with only cortical regions this can be only the
        cortical vector, other wise the whole data).
        :param permuted_imaging: imaging data permuted. Allows the user to
        specify the data to use (e.g., with only cortical regions this can be
        only the cortical vector, other wise the whole data).
        :param perm_indexes: indexes of the permuted data.
        :param scan_data: Original scam data, not zscored.
        :param gene_exp: gene expression data.
        :param gene_labels: gene labels.
        """

        def correlate(c1, c2):
            """Return the MATLAB style correlation between two vectors."""
            return np.corrcoef(np.hstack((c1, c2)), rowvar=False)[0, 1:]

        _res = pls_regression(gene_exp, imaging_data.reshape(
            imaging_data.shape[0], 1),
                              n_components=self.n_components,
                              n_boot=0, n_perm=0)
        r1 = correlate(_res.get("x_scores"), scan_data.reshape(
            scan_data.shape[0], 1))
        _weights = _res.get("x_weights")
        _scores = _res.get("x_scores")
        for i in range(r1.size):
            if r1[i] < 0:
                _weights[:, i] *= -1
                _scores[:, i] *= -1
        for _idx in range(self.n_components):
            _sort_weights_indexes = np.argsort(_weights[:, _idx],
                                               kind="mergesort")[::-1]
            self.orig.index[_idx, :] = _sort_weights_indexes
            self.orig.genes[_idx, :] = gene_labels[:, 0][_sort_weights_indexes]
            self.orig.weights[_idx, :] = _weights[:, _idx][
                _sort_weights_indexes]
            self.orig.zscored[_idx, :] = zscore(self.orig.weights[_idx, :],
                                                axis=0,
                                                ddof=1)
        for _iter in range(1000):
            _perm_imaging = permuted_imaging[:, _iter]
            _i_results = pls_regression(gene_exp, _perm_imaging.reshape(
                _perm_imaging.shape[0], 1),
                                        n_components=self.n_components,
                                        n_boot=0, n_perm=0)
            _weights_i = _i_results.get("x_weights")
            for _comp in range(self.n_components):
                _temp_weights = _weights_i[:, _comp]
                _new_weights = _temp_weights[self.orig.index[_comp, :]]
                _temp_genes = self.orig.weights[_comp, :]
                _corr = correlate(
                    _temp_genes.reshape(_temp_genes.size, 1),
                    _new_weights.reshape(_new_weights.shape[0], 1)
                )
                if _corr < 0:
                    _new_weights *= -1
                self.boot.weights[_comp, :, _iter] = _new_weights
        return

    def compute(self):
        """ Compute the p-values of the z-scored weights.
        The compute function calculates the zscores based on the original
        weights, orders the weights in descending order and then calculates
        the p-values and corrects for multiple comparisons using the
        Benjamini-Hochberg method.

        """
        for component in range(self.n_components):
            # Calculate the standard deviation of all the wights
            self.boot.std[component, :] = self.boot.weights[component, :,
                                          :].std(axis=1, ddof=1)
            # Calculate the zscore and store it sorted in descending order
            _z = self.orig.weights[component, :] / self.boot.std[component, :]
            self.boot.z_score[component, :] = np.sort(
                _z, axis=0, kind='mergesort')[::-1]
            _index = np.argsort(_z, axis=0, kind='mergesort')[::-1]
            # Reorder the genes according to the zscore
            self.boot.genes[component, :] = self.orig.genes[component, _index]
            # Calculate pvalue and pvalue corrected
            _p_val = norm.sf(abs(self.boot.z_score[component, :]))
            self.boot.pval[component, :] = _p_val
            _, _p_corr, _, _ = multipletests(_p_val, method='fdr_bh')
            self.boot.pval_corr[component, :] = _p_corr
        return

    def gsea(self, gene_set="lake", outdir=None):
        """Perform a GSEA analysis on the z-scored weights."""
        assert isinstance(self.orig, OrigPLS)
        assert isinstance(self.boot, BootPLS)
        gene_set = get_geneset(gene_set)
        for _component in range(self.n_components):
            gene_list = [gene for gene in self.orig.genes[
                                          _component, :].to_list()]
            rnk = pd.DataFrame(zip(gene_list,
                                   self.orig.z_score[_component, :]))
            gsea_results = gseapy.prerank(rnk, gene_set,
                                          outdir=None,
                                          seed=1234,
                                          permutation_num=100)
            _origin_es = gsea_results.res2d.es.to_numpy()
            _boot_es = np.zeros((_origin_es.shape[0], 1000))
            for i in range(1000):
                rnk = pd.DataFrame(zip(gene_list,
                                       zscore(
                                           self.boot.weights[_component, :, i],
                                           ddof=1)
                                       )
                                   )
                gsea_results = gseapy.prerank(rnk, gene_set,
                                              outdir=None,
                                              seed=1234,
                                              permutation_num=100)
                _boot_es[:, i] = gsea_results.res2d.es.to_numpy()
            _p_val = np.zeros((_origin_es.shape[0],))
            for i in range(_origin_es.shape[0]):
                _p_val[i] = np.sum(_boot_es[i, :] >= _origin_es[i]) / 1000
            # calculate the p-value corrected
            _, _p_corr, _, _ = multipletests(_p_val, method='fdr_bh')
            # Prepare data to save
            print("Prepare data...")
            _out_data = OrderedDict()
            _out_data["Term"] = gsea_results.res2d.axes[0].to_list()
            _out_data["es"] = gsea_results.res2d.values[:, 0]
            _out_data["nes"] = gsea_results.res2d.values[:, 1]
            _out_data["p_val"] = _p_val
            _out_data["fdr"] = _p_corr
            _out_data["genest_size"] = gsea_results.res2d.values[:, 4]
            _out_data["matched_size"] = gsea_results.res2d.values[:, 5]
            _out_data["matched_genes"] = gsea_results.res2d.values[:, 6]
            _out_data["ledge_genes"] = gsea_results.res2d.values[:, 7]
            out_df = pd.DataFrame.from_dict(_out_data)
            # TODO: make the output and clean the .csv file
            if outdir is not None:
                outdir = Path(outdir)
                assert outdir.exists()
                out_df.to_csv(
                    outdir / f"gsea_pls{_component + 1}_results.tsv",
                    index=False,
                    sep="\t")
                for _i in range(len(gsea_results.res2d.index)):
                    term = gsea_results.res2d.index[i]
                    gsea_results.results[term]["pval"] = _p_val[_i]
                    gsea_results.results[term]["fdr"] = _p_corr[_i]
                    gseaplot(rank_metric=gsea_results.ranking,
                             term=term,
                             **gsea_results.results[term],
                             ofname=f"{outdir}/imt_{term}_pls"
                                    f"{_component + 1}_prerank.pdf")


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
        self.genes = np.zeros((n_components, n_genes), dtype=object)
        self.index = np.zeros((n_components, n_genes), dtype=np.int32)
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
        self.genes = np.zeros((n_components, n_genes), dtype=object)
        self.std = np.zeros((n_components, n_genes))
        self._z_score = np.zeros((n_components, n_genes))
        self.pval = np.zeros((n_components, n_genes))
        self.pval_corr = np.zeros((n_components, n_genes))

    @property
    def z_score(self):
        """The z-scored weights.
        """
        return self._z_score


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
    def bootstrap_correlation(self, imaging_data, permuted_imaging,
                              gene_exp, gene_labels):
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
        :param list gene_labels: the labels of the genes.
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
                                            _ind, chunksize=15633)):
            self.gene_results.results.boot_corr.flat[ind] = res
        self.gene_results.results.genes = gene_labels
        self.gene_results.results.sort_genes()
        self.gene_results.results.compute_pval()
        return

    def gsea(self, gene_set="lake", outdir=None):
        """Perform GSEA on the correlation."""
        assert isinstance(self.gene_results.results, CorrGenes)
        print("Running GSEA...")
        gene_set = get_geneset(gene_set)
        # prepare the gene_list as a list of strings
        gene_list = [
            gene for gene in self.gene_results.results.genes[:, 0].tolist()
        ]
        # perform the GSEA on original results
        rnk = pd.DataFrame(
                zip(gene_list, self.gene_results.results.corr[0, :]))
        gsea_results = gseapy.prerank(rnk, gene_set,
                                      outdir=None,
                                      permutation_num=100,
                                      seed=1234)
        _origin_es = gsea_results.res2d.es.to_numpy()
        _boot_es = np.zeros((_origin_es.shape[0], 1000))
        # perform the GSEA on the permutations
        for i in range(1000):
            rnk = pd.DataFrame(
                zip(gene_list, self.gene_results.results.boot_corr[:, i])
            )
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
        # Prepare data to save
        print("Prepare data...")
        _out_data = OrderedDict()
        _out_data["Term"] = gsea_results.res2d.axes[0].to_list()
        _out_data["es"] = gsea_results.res2d.values[:, 0]
        _out_data["nes"] = gsea_results.res2d.values[:, 1]
        _out_data["p_val"] = _p_val
        _out_data["fdr"] = _p_corr
        _out_data["genest_size"] = gsea_results.res2d.values[:, 4]
        _out_data["matched_size"] = gsea_results.res2d.values[:, 5]
        _out_data["matched_genes"] = gsea_results.res2d.values[:, 6]
        _out_data["ledge_genes"] = gsea_results.res2d.values[:, 7]
        out_df = pd.DataFrame.from_dict(_out_data)
        # TODO: make the output and clean the .csv file
        if outdir is not None:
            outdir = Path(outdir)
            assert outdir.exists()
            out_df.to_csv(outdir / "gsea_results.tsv", index=False, sep="\t")
            for i in range(len(gsea_results.res2d.index)):
                term = gsea_results.res2d.index[i]
                gsea_results.results[term]["pval"] = _p_val[i]
                gsea_results.results[term]["fdr"] = _p_corr[i]
                gseaplot(rank_metric=gsea_results.ranking,
                         term=term,
                         **gsea_results.results[term],
                         ofname=f"{outdir}/imt_{term}_prerank.pdf")

    # --------- SAVE FUNCTIONS --------- #
    def save_results(self, outdir):
        """Save the results to a file."""
        outdir = Path(outdir)
        assert outdir.exists()
        # Create the data to save
        data_to_save = zip(self.gene_results.results.genes[:, 0],
                           self.gene_results.results.corr[0, :],
                           self.gene_results.results.pval[0, :],
                           self.gene_results.results.pval_corr[0, :])
        # Save the data
        pd.DataFrame(
            data_to_save, columns=["Gene", "Corr", "Pval", "Pval_corr"]
        ).to_csv(outdir / "corr_genes.tsv", index=False, sep='\t')
        return


def get_geneset(gene_set: str):
    """Returns the path to the geneset file, if it is not one of those
    defined in gseapy.

    :param str gene_set: the name of the geneset. If the geneset is "lake"
    or "pooled", the genesets are provided with the package data, otherwise
    the gene sets are from the gseapy package.
    """
    if gene_set.lower() == "lake":
        return str(Path(__file__).parent / "data" / "geneset_LAKE.gmt")
    elif gene_set.lower() == "pooled":
        return str(Path(__file__).parent / "data" /
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
        from the list of bootstrapped correlations.
        """
        # This calculation assumes that the order of the genes is the same
        # in both the original and the bootstrapped list. IF one is ordered,
        # make sure the order of the other is the same.
        print("compute pval")
        for i in range(self.n_genes):
            self.pval[0, i] = np.sum(
                self.boot_corr[i, :] >= self.corr[0, i]) / self._n_iter
        _, p_corr, _, _ = multipletests(self.pval[0, :], method='fdr_bh',
                                        is_sorted=False)
        self.pval_corr[0, :] = p_corr
        return

    @property
    def is_sorted(self):
        """Check if the list of genes is sorted.

        :return: True if the list of genes is sorted, False otherwise.
        """
        return False if self._index is None else True

    def sort_genes(self):
        """Order the genes in the list of genes. Both the order of the
        order of the bootstrapped genes are ordered.
        """
        print("sort genes")
        self._index = np.argsort(self.corr, axis=1, kind='mergesort')
        self.corr[0, :] = self.corr[0, self._index]
        self.genes[:] = self.genes[self._index]
        self.pval[0, :] = self.pval[0, self._index]
        self.pval_corr[0, :] = self.pval_corr[0, self._index]
        for i in range(self._n_iter):
            self.boot_corr[:, i] = self.boot_corr[:, i][self._index]
        return
