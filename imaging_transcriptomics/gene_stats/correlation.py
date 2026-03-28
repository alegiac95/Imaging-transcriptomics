from __future__ import annotations

import warnings

import numpy as np

from .._logging import get_logger
from ..stats_utils import bh_fdr, empirical_signed_pvalues, max_t_fwer_abs, minimum_bh_resolution

logger = get_logger("genes")


class CorrGenes:
    """Store gene-wise statistics for the correlation workflow."""

    def __init__(self, n_iter=1000, n_genes=None, store_boot_corr: bool = True):
        """Create storage for observed and permuted correlation statistics."""

        self.n_genes = int(n_genes) if n_genes is not None else 15633
        self._n_iter = n_iter
        self.boot_corr = np.zeros((self.n_genes, self._n_iter), dtype=np.float32) if store_boot_corr else None
        self.corr = np.zeros((1, self.n_genes))
        self.genes = np.zeros((self.n_genes, 1), dtype=object)
        self.pval = np.zeros((1, self.n_genes))
        self.pval_corr = np.zeros((1, self.n_genes))
        self.pval_fwer = np.zeros((1, self.n_genes))
        self._index = None
        self._signed_counts = np.zeros(self.n_genes, dtype=np.int64)
        self._max_t_counts = np.zeros(self.n_genes, dtype=np.int64)
        self._processed_iterations = 0

    def accumulate_boot_corr(self, boot_corr: np.ndarray, *, start: int | None = None) -> None:
        """Accumulate one chunk of permuted correlation nulls."""

        chunk = np.asarray(boot_corr, dtype=float)
        if chunk.ndim != 2 or chunk.shape[0] != self.n_genes:
            raise ValueError("boot_corr must be a 2D matrix with one row per gene.")
        if self.boot_corr is not None:
            if start is None:
                raise ValueError("start is required when storing boot_corr chunks.")
            self.boot_corr[:, start : start + chunk.shape[1]] = chunk
        observed = self.corr[0, :]
        pos_counts = np.sum(chunk >= observed.reshape(-1, 1), axis=1)
        neg_counts = np.sum(chunk <= observed.reshape(-1, 1), axis=1)
        self._signed_counts += np.where(observed >= 0, pos_counts, neg_counts)
        perm_max_abs = np.max(np.abs(chunk), axis=0)
        self._max_t_counts += np.sum(
            perm_max_abs.reshape(1, -1) >= np.abs(observed).reshape(-1, 1),
            axis=1,
        )
        self._processed_iterations += chunk.shape[1]

    def compute_pval(self):
        """Compute gene-wise nominal, BH-corrected, and maxT-corrected p-values."""

        logger.info("Computing p values.")
        if self._processed_iterations > 0:
            denominator = self._processed_iterations + 1
            self.pval[0, :] = (self._signed_counts + 1) / denominator
            self.pval_fwer[0, :] = (self._max_t_counts + 1) / denominator
            n_permutations = self._processed_iterations
        elif self.boot_corr is not None:
            self.pval[0, :] = empirical_signed_pvalues(self.corr[0, :], self.boot_corr)
            self.pval_fwer[0, :] = max_t_fwer_abs(self.corr[0, :], self.boot_corr)
            n_permutations = self._n_iter
        else:
            raise RuntimeError("No correlation nulls are available to compute p-values.")

        min_possible_p = 1.0 / (n_permutations + 1)
        min_possible_bh = minimum_bh_resolution(self.n_genes, n_permutations)
        if min_possible_bh >= 1.0:
            warnings.warn(
                "Correlation gene FDR uses Benjamini-Hochberg on permutation p-values, "
                f"but with {n_permutations} permutations across {self.n_genes} genes the "
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
        self._index = np.argsort(self.corr[0, :], kind="mergesort")[::-1]
        self.corr[0, :] = self.corr[0, self._index]
        self.genes = self.genes[self._index, :]
        self.pval[0, :] = self.pval[0, self._index]
        self.pval_corr[0, :] = self.pval_corr[0, self._index]
        self.pval_fwer[0, :] = self.pval_fwer[0, self._index]
        if self.boot_corr is not None:
            self.boot_corr = self.boot_corr[self._index, :]
        self._signed_counts = self._signed_counts[self._index]
        self._max_t_counts = self._max_t_counts[self._index]
        return
