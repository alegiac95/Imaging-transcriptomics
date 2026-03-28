from __future__ import annotations

from .correlation import CorrGenes
from .pls import PLSGenes


class GeneResults:
    """Expose a uniform view over correlation and PLS gene result containers."""

    def __init__(self, method, **kwargs):
        """Initialize the workflow-specific gene result container."""

        self.method = method
        n_genes = kwargs.get("n_genes")
        if self.method == "pls":
            self.results = PLSGenes(
                kwargs.get("n_components"),
                n_iter=kwargs.get("n_iter", 1000),
                n_genes=n_genes,
                store_weights=kwargs.get("store_weights", True),
            )
        elif self.method == "corr":
            self.results = CorrGenes(
                n_iter=kwargs.get("n_iter"),
                n_genes=n_genes,
                store_boot_corr=kwargs.get("store_boot_corr", True),
            )
        else:
            raise ValueError(f"The method {method} is not supported.")

    @property
    def n_genes(self):
        return self.results.n_genes

    @property
    def genes(self):
        if isinstance(self.results, PLSGenes):
            return self.results.orig.genes
        if isinstance(self.results, CorrGenes):
            return self.results.genes
        return None

    @property
    def scores(self):
        if isinstance(self.results, PLSGenes):
            return self.results.orig.weights
        if isinstance(self.results, CorrGenes):
            return self.results.corr
        return None

    @property
    def boot(self):
        if isinstance(self.results, PLSGenes):
            return self.results.boot.weights
        if isinstance(self.results, CorrGenes):
            return self.results.boot_corr
        return None

    @property
    def pvals(self):
        if isinstance(self.results, PLSGenes):
            return self.results.boot.pval
        if isinstance(self.results, CorrGenes):
            return self.results.pval
        return None

    @property
    def pvals_corr(self):
        if isinstance(self.results, PLSGenes):
            return self.results.boot.pval_corr
        if isinstance(self.results, CorrGenes):
            return self.results.pval_corr
        return None

    @property
    def pvals_fwer(self):
        if isinstance(self.results, PLSGenes):
            return self.results.boot.pval_fwer
        if isinstance(self.results, CorrGenes):
            return self.results.pval_fwer
        return None
