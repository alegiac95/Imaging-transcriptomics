import warnings
from pathlib import Path
from scipy.stats import zscore
import numpy as np
from collections import OrderedDict
import pandas as pd
from ._logging import get_logger
from .gsea_utils import (
    gsea_style_fdr,
    make_prerank_table,
    nominal_pvalues_from_nulls,
    normalize_enrichment_nulls,
    normalize_enrichment_scores,
    run_prerank,
)
from .genesets import get_geneset
from .ora import ora_from_gene_table
from .pls_backend import pls_regression
from .stats_utils import bh_fdr, empirical_signed_pvalues, max_t_fwer_abs, minimum_bh_resolution, two_sided_z_pvalues

logger = get_logger(__name__)


def _correlate_pls_scores(scores: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Correlate each PLS score vector with the input regional imaging values."""

    stacked = np.hstack((np.asarray(scores, dtype=float), np.asarray(values, dtype=float).reshape(-1, 1)))
    return np.corrcoef(stacked, rowvar=False)[0, 1:]


def _rowwise_corrsign(
    reference: np.ndarray,
    candidate: np.ndarray,
    reference_centered: np.ndarray,
    reference_ss: np.ndarray,
) -> np.ndarray:
    """Return sign flips that align candidate rows to reference rows by correlation."""

    centered = candidate - candidate.mean(axis=1, keepdims=True)
    denom = np.sqrt(reference_ss * np.sum(centered * centered, axis=1))
    corr = np.divide(
        np.sum(reference_centered * centered, axis=1),
        denom,
        out=np.zeros(reference.shape[0], dtype=float),
        where=denom != 0,
    )
    return np.where(corr < 0, -1.0, 1.0)


# --------- GENE ANALYSIS --------- #
class GeneResults:
    """Expose a uniform view over correlation and PLS gene result containers."""

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
        n_genes = kwargs.get("n_genes")
        if self.method == "pls":
            self.results = PLSGenes(kwargs.get("n_components"),
                                    n_iter=kwargs.get("n_iter", 1000),
                                    n_genes=n_genes)
        elif self.method == "corr":
            self.results = CorrGenes(n_iter=kwargs.get("n_iter"), n_genes=n_genes)
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

    @property
    def pvals_fwer(self):
        if isinstance(self.results, PLSGenes):
            return self.results.boot.pval_fwer
        if isinstance(self.results, CorrGenes):
            return self.results.pval_fwer
        return None


# --------- PLS GENES --------- #
class PLSGenes:
    """Store original and permutation-based gene statistics for PLS results."""

    def __init__(self, n_components, n_iter=1000, n_genes=None):
        """ Initialize the results of the PLS analysis. The result will
        include both the permuted and the original results. The class
        contains two subclasses, one for the original results and one for
        the bootrapped results.

        :param int n_components: number of components used for the analysis.
        """
        self.n_genes = int(n_genes) if n_genes is not None else 15633
        self.n_components = n_components
        self.n_iter = n_iter
        self.orig = OrigPLS(n_components, self.n_genes)
        self.boot = BootPLS(n_components, self.n_genes, n_iter=n_iter)
        self._orig_centered = None
        self._orig_ss = None

    def prepare_from_fit(self, fit_result, scan_data, gene_labels):
        """Align, sort, and z-score original PLS gene weights per component."""

        weights = np.asarray(fit_result.get("x_weights"), dtype=float).copy()
        scores = np.asarray(fit_result.get("x_scores"), dtype=float).copy()
        score_corr = _correlate_pls_scores(scores, scan_data)
        for component in range(score_corr.size):
            if score_corr[component] < 0:
                weights[:, component] *= -1
                scores[:, component] *= -1
        for component in range(self.n_components):
            sort_index = np.argsort(weights[:, component], kind="mergesort")[::-1]
            self.orig.index[component, :] = sort_index
            self.orig.genes[component, :] = gene_labels[:, 0][sort_index]
            self.orig.weights[component, :] = weights[:, component][sort_index]
            self.orig.zscored[component, :] = zscore(self.orig.weights[component, :], axis=0, ddof=1)
        self._orig_centered = self.orig.weights - self.orig.weights.mean(axis=1, keepdims=True)
        self._orig_ss = np.sum(self._orig_centered * self._orig_centered, axis=1)
        return

    def store_permuted_weights(self, iteration: int, x_weights):
        """Store one permuted PLS fit after aligning it to the original gene order."""

        if self._orig_centered is None or self._orig_ss is None:
            raise RuntimeError("Call prepare_from_fit() before storing permuted weights.")
        weights = np.asarray(x_weights, dtype=float).T
        reordered = np.take_along_axis(weights, self.orig.index, axis=1)
        sign = _rowwise_corrsign(self.orig.weights, reordered, self._orig_centered, self._orig_ss)
        self.boot.weights[:, :, iteration] = reordered * sign.reshape(-1, 1)
        return

    def boot_genes(self, imaging_data, permuted_imaging,
                   scan_data, gene_exp, gene_labels):
        """Legacy helper to refit PLS for every permutation and store gene weights."""
        logger.info("Performing bootstrapping of the genes.")

        _res = pls_regression(gene_exp, imaging_data.reshape(
            imaging_data.shape[0], 1),
                              n_components=self.n_components,
                              n_boot=0, n_perm=0)
        self.prepare_from_fit(_res, scan_data, gene_labels)
        if permuted_imaging.shape[1] != self.boot.weights.shape[2]:
            raise ValueError("The number of bootstrapped permutations does "
                             "not match the configured iteration count.")
        for _iter in range(self.boot.weights.shape[2]):
            _perm_imaging = permuted_imaging[:, _iter]
            _i_results = pls_regression(gene_exp, _perm_imaging.reshape(
                                        _perm_imaging.shape[0], 1),
                                        n_components=self.n_components,
                                        n_boot=0, n_perm=0)
            self.store_permuted_weights(_iter, _i_results.get("x_weights"))
        return

    def compute(self):
        """Compute sorted PLS gene statistics from original and permuted weights.

        This step estimates a bootstrap standard deviation per gene and
        component, derives z-scores from the original weights, computes
        two-sided z-based nominal p-values, then adds BH FDR and maxT-style
        FWER correction before sorting the tables in descending weight order.
        """
        logger.info("Calculating statistics.")
        self.boot.std[:, :] = self.boot.weights.std(axis=2, ddof=1)
        safe_std = np.where(self.boot.std == 0, np.finfo(float).eps, self.boot.std)
        zscores = self.orig.weights / safe_std
        indices = np.argsort(zscores, axis=1, kind='mergesort')[:, ::-1]
        raw_pval = two_sided_z_pvalues(zscores)
        raw_pval_fwer = np.zeros((self.n_components, self.n_genes), dtype=float)
        self.boot.weights_sorted[:, :] = np.take_along_axis(self.orig.weights, indices, axis=1)
        self.boot.z_score[:, :] = np.take_along_axis(zscores, indices, axis=1)
        self.boot.genes[:, :] = np.take_along_axis(self.orig.genes, indices, axis=1)
        for component in range(self.n_components):
            raw_pval_fwer[component, :] = max_t_fwer_abs(
                self.orig.weights[component, :],
                self.boot.weights[component, :, :],
            )
            corrected = bh_fdr(raw_pval[component, :])
            self.boot.pval[component, :] = raw_pval[component, indices[component, :]]
            self.boot.pval_corr[component, :] = corrected[indices[component, :]]
            self.boot.pval_fwer[component, :] = raw_pval_fwer[component, indices[component, :]]
        return

    def gsea(self, gene_set="lake", outdir=None, gene_limit=1500, n_iter=1000):
        """Run preranked GSEA on the PLS gene ranking for each component.

        Observed ES values come from the original component z-scores, while NES,
        nominal p-values, and q-values are recalculated from the external null
        built from permuted component weights.
        """
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
            rnk = make_prerank_table(
                gene_list,
                self.orig.zscored[_component, :],
            )
            gsea_results = run_prerank(gseapy, rnk, gene_set,
                                       max_size=gene_limit,
                                       outdir=None,
                                       seed=1234,
                                       permutation_num=0)
            _origin_es = gsea_results.res2d.es.to_numpy()
            _boot_es = np.zeros((_origin_es.shape[0], n_iter))
            for i in range(n_iter):
                rnk = make_prerank_table(
                    gene_list,
                    zscore(
                        self.boot.weights[_component, :, i],
                        ddof=1,
                    ),
                )
                gsea_res = run_prerank(gseapy, rnk, gene_set,
                                       max_size=gene_limit,
                                       outdir=None,
                                       seed=1234,
                                       permutation_num=0)
                _boot_es[:, i] = gsea_res.res2d.es.to_numpy()
            _nes = normalize_enrichment_scores(_origin_es, _boot_es)
            _nes_null = normalize_enrichment_nulls(_origin_es, _boot_es)
            _p_val = nominal_pvalues_from_nulls(_origin_es, _boot_es)
            _p_corr = gsea_style_fdr(_nes, _nes_null)
            # Prepare data to save
            _out_data = OrderedDict()
            _out_data["Term"] = gsea_results.res2d.axes[0].to_list()
            _out_data["es"] = gsea_results.res2d.values[:, 0]
            _out_data["nes"] = _nes
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

    def ora(self, gene_set="lake", outdir=None, p_threshold=0.05):
        """Run ORA on positive and negative component gene tails separately."""

        assert isinstance(self.orig, OrigPLS)
        assert isinstance(self.boot, BootPLS)
        logger.info("Performing ORA.")
        results: list[dict[str, pd.DataFrame]] = []
        for _component in range(self.n_components):
            gene_table = pd.DataFrame(
                {
                    "gene": self.boot.genes[_component, :],
                    "zscore": self.boot.z_score[_component, :],
                    "p_value": self.boot.pval[_component, :],
                }
            )
            ora_tables = ora_from_gene_table(
                gene_table,
                gene_set=gene_set,
                score_column="zscore",
                p_threshold=p_threshold,
            )
            if outdir is not None:
                logger.info("Saving ORA results for PLS component %d.", _component + 1)
                output_dir = Path(outdir)
                assert output_dir.exists()
                for direction, table in ora_tables.items():
                    table.to_csv(
                        output_dir / f"ora_pls{_component + 1}_{direction}.tsv",
                        index=False,
                        sep="\t",
                    )
            results.append(ora_tables)
        return results


# --------- ORIG PLS --------- #
class OrigPLS:
    """Hold original per-component PLS gene rankings and summary statistics."""

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
    """Hold permutation-derived PLS gene weights and correction outputs."""

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
        self.weights_sorted = np.zeros((n_components, n_genes))
        self.std = np.zeros((n_components, n_genes))
        self._z_score = np.zeros((n_components, n_genes))
        self.pval = np.zeros((n_components, n_genes))
        self.pval_corr = np.zeros((n_components, n_genes))
        self.pval_fwer = np.zeros((n_components, n_genes))

    @property
    def z_score(self):
        """Return the sorted per-gene z-scores for each PLS component."""
        return self._z_score


# --------- CORRELATION  GENES  --------- #
class CorrGenes:
    """Store gene-wise statistics for the correlation workflow."""
    def __init__(self, n_iter=1000, n_genes=None):
        """Create storage for observed and permuted correlation statistics."""
        self.n_genes = int(n_genes) if n_genes is not None else 15633
        self._n_iter = n_iter
        self.boot_corr = np.zeros((self.n_genes, self._n_iter))
        self.corr = np.zeros((1, self.n_genes))
        self.genes = np.zeros((self.n_genes, 1), dtype=object)
        self.pval = np.zeros((1, self.n_genes))
        self.pval_corr = np.zeros((1, self.n_genes))
        self.pval_fwer = np.zeros((1, self.n_genes))
        self._index = None

    def compute_pval(self):
        """Compute gene-wise nominal, BH-corrected, and maxT-corrected p-values."""
        # This calculation assumes that the order of the genes is the same
        # in both the original and the bootstrapped list. IF one is ordered,
        # make sure the order of the other is the same.
        logger.info("Computing p values.")
        self.pval[0, :] = empirical_signed_pvalues(self.corr[0, :], self.boot_corr)
        self.pval_fwer[0, :] = max_t_fwer_abs(self.corr[0, :], self.boot_corr)
        min_possible_p = 1.0 / (self._n_iter + 1)
        min_possible_bh = minimum_bh_resolution(self.n_genes, self._n_iter)
        if min_possible_bh >= 1.0:
            warnings.warn(
                "Correlation gene FDR uses Benjamini-Hochberg on permutation p-values, "
                f"but with {self._n_iter} permutations across {self.n_genes} genes the "
                f"smallest possible nominal p-value is {min_possible_p:.6g}, so adjusted "
                "p-values will collapse to 1. Increase permutations above the number of "
                "genes for non-trivial gene-level FDR.",
                RuntimeWarning,
                stacklevel=2,
            )
        self.pval_corr[0, :] = bh_fdr(self.pval[0, :])
        return

    @property
    def is_sorted(self):
        """Whether the gene table has already been sorted by observed score."""
        return self._index is not None

    def sort_genes(self):
        """Sort observed and permuted gene statistics by descending correlation."""
        logger.info("Sorting genes in descending order.")
        self._index = np.argsort(self.corr[0, :], kind='mergesort')[::-1]
        self.corr[0, :] = self.corr[0, self._index]
        self.genes = self.genes[self._index, :]
        self.pval[0, :] = self.pval[0, self._index]
        self.pval_corr[0, :] = self.pval_corr[0, self._index]
        self.pval_fwer[0, :] = self.pval_fwer[0, self._index]
        self.boot_corr = self.boot_corr[self._index, :]
        return
