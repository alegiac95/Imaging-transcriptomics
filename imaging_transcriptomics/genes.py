import logging
import logging.config
import yaml
from pathlib import Path

from scipy.stats import zscore, norm
from statsmodels.stats.multitest import multipletests
import numpy as np


cfg_file_path = Path(__file__).parent / "log_config.yaml"
with open(cfg_file_path, "r") as config_file:
    log_cfg = yaml.safe_load(config_file.read())

logging.config.dictConfig(log_cfg)
logger = logging.getLogger("genes")
logger.setLevel(logging.DEBUG)


class OriginalResults(dict):
    """Class to hold the result from gene analysis (non bootstrapped)."""

    def __init__(self, n_comp):
        super().__init__()
        self.pls_weights = [None] * n_comp
        self.pls_gene = [None] * n_comp
        self.gene_id = [None] * n_comp
        self.index = [None] * n_comp
        self.pls_weights_z = [None] * n_comp
        logger.debug("OriginalResults class initialized with %s components", n_comp)

    def set_result_values(self, idx, pls_weights, x, pls_genes, gene_id):
        """Set the values for all class attributes.

        :param idx: index of the component to set the values in the array (must be in the range 0:n_components-1).
        :param pls_weights: weights to pass to the array at the index idx.
        :param pls_genes: genes to pass to the array at the index idx
        :param gene_id: id of the genes to pass to the array at the index idx
        :param x: index of the pls genes in the original gene list.
        :return:
        """
        self.pls_weights[idx - 1] = np.array(pls_weights)
        logger.debug("Weights at index %s set as %s", idx - 1, pls_weights)
        self.pls_gene[idx - 1] = np.array(pls_genes)
        logger.debug("Genes at index %s set as %s", idx - 1, pls_genes)
        self.gene_id[idx - 1] = gene_id
        logger.debug("Genes ids at index %s set as %s", idx - 1, gene_id)
        self.index[idx - 1] = np.array(x)
        self.pls_weights_z[idx - 1] = np.array(zscore(pls_weights, axis=0, ddof=1))


class BootResults(dict):
    """Class to hold the results from the bootstrapping gene analysis."""

    def __init__(self, n_comp, dim1, dim2, n_perm=1000):
        super().__init__()
        self.pls_weights_boot = [np.zeros((dim1, dim2, n_perm))] * n_comp
        self.std = [None] * n_comp
        self.z_scores = [None] * n_comp
        self.pls_genes = [None] * n_comp
        self.pval = [None] * n_comp
        self.pval_corrected = [None] * n_comp
        logger.debug("BootResult class initialised with %s components", n_comp)

    def compute_values(self, n_comp, original_weights, original_ids):
        """Compute the values of the bootstrap for each of the components.

        :param n_comp: number of PLS components.
        :param original_weights: weights obtained from the original analysis (not bootstrapped)
        :param original_ids: original ids (labels) from the original analysis (not bootstrapped)
        :return:
        """
        logger.info("Computing bootstrap gene results.")
        for component in range(1, n_comp + 1):
            logger.debug("Computing results for component %s", component)
            self.std[component - 1] = self.pls_weights_boot[component - 1][
                :, component - 1, :
            ].std(ddof=1, axis=1)
            __temp = original_weights[component - 1] / self.std[component - 1]
            self.z_scores[component - 1] = np.sort(__temp, kind="mergesort")[::-1]
            __idx = np.argsort(__temp, kind="mergesort")[::-1]
            self.pls_genes[component - 1] = original_ids[component - 1][__idx]
            __p = norm.sf(abs(self.z_scores[component - 1]))
            self.pval[component - 1] = __p
            _, __p_corr, _, _ = multipletests(
                __p[::-1].reshape(1, 15633), method="fdr_bh", is_sorted=True
            )
            self.pval_corrected[component - 1] = __p_corr

    def compute_correlation(self, original_corr, originals_ids, original_index):
        """Compute the p value of the correlation after the permutation."""
        logger.info("Computing bootstrap correlation results")
        n_iter = 1000
        for i in range(n_iter):
            tmp = self.pls_weights_boot[:, i]
            self.pls_weights_boot[:, i] = tmp[original_index]
        p = np.zeros(15633)
        for i in range(15633):
            original = original_corr[i]
            boot = self.pls_weights_boot[i, :]
            p[i] = float(len(boot[np.nonzero(boot >= original)])) / n_iter
        self.pval = p
        _, p_corrected, _, _ = multipletests(p, method="fdr_bh",
                                             is_sorted=False)
        self.pval_corrected = p_corrected
        self.pls_genes = np.array(originals_ids)


class GeneResults(dict):
    """Class to save the gene results."""

    def __init__(self, n_comp, dim1, dim2):
        super().__init__()
        self.original_results = OriginalResults(n_comp)
        self.boot_results = BootResults(n_comp, dim1, dim2)
