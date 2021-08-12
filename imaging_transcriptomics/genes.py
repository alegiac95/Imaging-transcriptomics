from scipy.stats import zscore, norm
from statsmodels.stats.multitest import multipletests
import numpy as np


class OriginalResults(dict):
    def __init__(self, n_comp):
        super().__init__()
        self.pls_weights = [None] * n_comp
        self.pls_gene = [None] * n_comp
        self.gene_id = [None] * n_comp
        self.index = [None] * n_comp
        self.pls_weights_z = [None] * n_comp

    def set_result_values(self, idx, pls_weights, x, pls_genes, gene_id):
        """Set the values for all class attributes.
        
        :param idx: 
        :param pls_weights: 
        :param pls_genes: 
        :param gene_id: 
        :param x: 
        :return: 
        """
        self.pls_weights[idx - 1] = np.array(pls_weights),
        self.pls_gene[idx - 1] = np.array(pls_genes),
        self.gene_id[idx - 1] = gene_id,
        self.index[idx - 1] = np.array(x)
        self.pls_weights_z[idx - 1] = np.array(zscore(pls_weights, axis=0, ddof=1))


class BootResults(dict):
    def __init__(self, n_comp, dim1, dim2, n_perm=1000):
        super().__init__()
        self.pls_weights_boot = [np.zeros((dim1, dim2, n_perm))] * n_comp
        self.std = [None] * n_comp
        self.z_scores = [None] * n_comp
        self.pls_genes = [None] * n_comp
        self.pval = [None] * n_comp
        self.pval_corrected = [None] * n_comp

    def compute_values(self, n_comp, original_weights, original_ids):
        """Compute the values of the bootstrap for each of the components.

        :param n_comp:
        :param original_weights:
        :param original_ids:
        :return:
        """
        for component in range(1, n_comp + 1):
            self.std[component - 1] = self.pls_weights_boot[component - 1][:, component - 1, :].std(ddof=1, axis=1)
            __temp = original_weights[component - 1] / self.std[component - 1]
            self.z_scores[component - 1] = np.sort(__temp, kind='mergesort')[::-1]
            __idx = np.argsort(__temp, kind='mergesort')[::-1]
            self.pls_genes[component - 1] = original_ids[component - 1][0][__idx]
            __p = norm.sf(abs(self.z_scores[component - 1]))
            self.pval[component - 1] = __p
            _, __p_corr, _, _ = multipletests(__p[::-1].reshape(1, 15633),
                                              method="fdr_bh", is_sorted=True)
            self.pval_corrected[component - 1] = __p_corr


class GeneResults(dict):
    """Class to save the gene results."""
    def __init__(self, n_comp, dim1, dim2):
        super().__init__()
        self.original_results = OriginalResults(n_comp)
        self.boot_results = BootResults(n_comp, dim1, dim2)
