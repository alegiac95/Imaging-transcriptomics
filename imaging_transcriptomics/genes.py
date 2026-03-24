import logging
from pathlib import Path
from scipy.stats import zscore, norm
from statsmodels.stats.multitest import multipletests
import numpy as np
from collections import OrderedDict
import pandas as pd
from ._logging import get_logger
from .genesets import get_geneset
from .pls_backend import pls_regression

logger = get_logger("genes")
logger.setLevel(logging.DEBUG)


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
            self.results = PLSGenes(kwargs.get("n_components"),
                                    n_iter=kwargs.get("n_iter", 1000))
        elif self.method == "corr":
            self.results = CorrGenes(n_iter=kwargs.get("n_iter"))
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


# --------- PLS GENES --------- #
class PLSGenes:
    def __init__(self, n_components, n_iter=1000):
        """ Initialize the results of the PLS analysis. The result will
        include both the permuted and the original results. The class
        contains two subclasses, one for the original results and one for
        the bootrapped results.

        :param int n_components: number of components used for the analysis.
        """
        self.n_genes = 15633
        self.n_components = n_components
        self.n_iter = n_iter
        self.orig = OrigPLS(n_components, self.n_genes)
        self.boot = BootPLS(n_components, self.n_genes, n_iter=n_iter)

    def boot_genes(self, imaging_data, permuted_imaging,
                   scan_data, gene_exp, gene_labels):
        """Bootstrapping on the PLS components.

        :param imaging_data: imaging data. Allows the user to specify the
        data to use (e.g., with only cortical regions this can be only the
        cortical vector, other wise the whole data).
        :param permuted_imaging: imaging data permuted. Allows the user to
        specify the data to use (e.g., with only cortical regions this can be
        only the cortical vector, other wise the whole data).
        :param scan_data: Original scam data, not zscored.
        :param gene_exp: gene expression data.
        :param gene_labels: gene labels.
        """
        logger.info("Performing bootstrapping of the genes.")

        def correlate(c1, c2):
            """Return the MATLAB style correlation between two vectors."""
            return np.corrcoef(np.hstack((c1, c2)), rowvar=False)[0, 1:]

        def rowwise_corrsign(reference, candidate, reference_centered, reference_ss):
            centered = candidate - candidate.mean(axis=1, keepdims=True)
            denom = np.sqrt(reference_ss * np.sum(centered * centered, axis=1))
            corr = np.divide(
                np.sum(reference_centered * centered, axis=1),
                denom,
                out=np.zeros(reference.shape[0], dtype=float),
                where=denom != 0,
            )
            return np.where(corr < 0, -1.0, 1.0)

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
        orig_centered = self.orig.weights - self.orig.weights.mean(axis=1, keepdims=True)
        orig_ss = np.sum(orig_centered * orig_centered, axis=1)
        if permuted_imaging.shape[1] != self.boot.weights.shape[2]:
            raise ValueError("The number of bootstrapped permutations does "
                             "not match the configured iteration count.")
        for _iter in range(self.boot.weights.shape[2]):
            _perm_imaging = permuted_imaging[:, _iter]
            _i_results = pls_regression(gene_exp, _perm_imaging.reshape(
                                        _perm_imaging.shape[0], 1),
                                        n_components=self.n_components,
                                        n_boot=0, n_perm=0)
            _weights_i = _i_results.get("x_weights").T
            _new_weights = np.take_along_axis(_weights_i, self.orig.index, axis=1)
            _sign = rowwise_corrsign(self.orig.weights, _new_weights, orig_centered, orig_ss)
            self.boot.weights[:, :, _iter] = _new_weights * _sign.reshape(-1, 1)
        return

    def compute(self):
        """ Compute the p-values of the z-scored weights.
        The compute function calculates the zscores based on the original
        weights, orders the weights in descending order and then calculates
        the p-values and corrects for multiple comparisons using the
        Benjamini-Hochberg method.
        """
        logger.info("Calculating statistics.")
        self.boot.std[:, :] = self.boot.weights.std(axis=2, ddof=1)
        safe_std = np.where(self.boot.std == 0, np.finfo(float).eps, self.boot.std)
        zscores = self.orig.weights / safe_std
        indices = np.argsort(zscores, axis=1, kind='mergesort')[:, ::-1]
        self.boot.z_score[:, :] = np.take_along_axis(zscores, indices, axis=1)
        self.boot.genes[:, :] = np.take_along_axis(self.orig.genes, indices, axis=1)
        self.boot.pval[:, :] = 2 * norm.sf(np.abs(self.boot.z_score))
        for component in range(self.n_components):
            _, self.boot.pval_corr[component, :], _, _ = multipletests(
                self.boot.pval[component, :],
                method='fdr_bh',
                is_sorted=False,
            )
        return

    def gsea(self, gene_set="lake", outdir=None, gene_limit=1500, n_iter=1000):
        """Perform a GSEA analysis on the z-scored weights."""
        assert isinstance(self.orig, OrigPLS)
        assert isinstance(self.boot, BootPLS)
        logger.info("Performing GSEA.")
        try:
            import gseapy
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError("gseapy is required to run GSEA analyses.") from exc
        if Path(gene_set).exists() and Path(gene_set).is_file() and Path(
                gene_set).suffix == ".gmt":
            gene_set = Path(gene_set)
        else:
            gene_set = get_geneset(gene_set)
        for _component in range(self.n_components):
            gene_list = list(self.orig.genes[_component, :])
            rnk = pd.DataFrame(zip(gene_list,
                                   self.orig.zscored[_component, :]))
            gsea_results = gseapy.prerank(rnk, gene_set,
                                          max_size=gene_limit,
                                          outdir=None,
                                          seed=1234,
                                          permutation_num=n_iter)
            _origin_es = gsea_results.res2d.es.to_numpy()
            _boot_es = np.zeros((_origin_es.shape[0], n_iter))
            for i in range(n_iter):
                rnk = pd.DataFrame(zip(gene_list,
                                       zscore(
                                           self.boot.weights[_component, :, i],
                                           ddof=1)
                                       )
                                   )
                gsea_res = gseapy.prerank(rnk, gene_set,
                                          max_size=gene_limit,
                                          outdir=None,
                                          seed=1234,
                                          permutation_num=1)
                _boot_es[:, i] = gsea_res.res2d.es.to_numpy()
            _p_val = np.zeros((_origin_es.shape[0],))
            for i in range(_origin_es.shape[0]):
                _count = np.sum(_boot_es[i, :] >= _origin_es[i]) \
                    if _origin_es[i] >= 0 else \
                    np.sum(_boot_es[i, :] <= _origin_es[i])
                _p_val[i] = (_count + 1) / (n_iter + 1)
            # calculate the p-value corrected
            _, _p_corr, _, _ = multipletests(_p_val, method='fdr_bh',
                               is_sorted=False)
            # Prepare data to save
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
            if outdir is not None:
                logger.info("Saving GSEA results.")
                outdir = Path(outdir)
                assert outdir.exists()
                out_df.to_csv(
                    outdir / f"gsea_pls{_component + 1}_results.tsv",
                    index=False,
                    sep="\t")


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
        self.n_iter = n_iter
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
        self.genes = np.zeros((self.n_genes, 1), dtype=object)
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
        logger.info("Computing p values.")
        pos_counts = np.sum(self.boot_corr >= self.corr.T, axis=1)
        neg_counts = np.sum(self.boot_corr <= self.corr.T, axis=1)
        counts = np.where(self.corr[0, :] >= 0, pos_counts, neg_counts)
        self.pval[0, :] = (counts + 1) / (self._n_iter + 1)
        _, p_corr, _, _ = multipletests(self.pval[0, :], method='fdr_bh',
                                        is_sorted=False)

        self.pval_corr[0, :] = p_corr
        return

    @property
    def is_sorted(self):
        """Check if the list of genes is sorted.

        :return: True if the list of genes is sorted, False otherwise.
        """
        return self._index is not None

    def sort_genes(self):
        """Order the genes in the list of genes. Both the order of the
        order of the bootstrapped genes are ordered.
        """
        logger.info("Sorting genes in descending order.")
        self._index = np.argsort(self.corr[0, :], kind='mergesort')[::-1]
        self.corr[0, :] = self.corr[0, self._index]
        self.genes = self.genes[self._index, :]
        self.pval[0, :] = self.pval[0, self._index]
        self.pval_corr[0, :] = self.pval_corr[0, self._index]
        self.boot_corr = self.boot_corr[self._index, :]
        return
