from __future__ import annotations

import logging

import numpy as np

from ._logging import get_logger
from .genes import GeneResults
from .pls_backend import pls_regression

logger = get_logger("genes")
logger.setLevel(logging.DEBUG)


class PLSAnalysis:
    """Perform PLS regression for the imaging transcriptomics workflows."""

    def __init__(self, imaging_data, gene_exp, n_components: int, var: float, n_iter: int = 1000):
        self.var, self.n_components, self.components_var = self.set_coef(
            imaging_data,
            gene_exp,
            n_components=n_components,
            var=var,
        )
        self.n_iter = n_iter
        self._p_val = np.zeros(self.n_components, dtype=float)
        self._r2 = np.zeros(self.n_components, dtype=float)
        self.gene_results = GeneResults(
            "pls",
            n_components=self.n_components,
            n_iter=n_iter,
            n_genes=np.asarray(gene_exp).shape[1],
        )

    @staticmethod
    def check_var(var: float):
        if var < 0:
            raise ValueError("The variance must be a positive number.")
        if var > 1:
            raise ValueError("The variance must be a number between 0 and 1.")
        return var

    @staticmethod
    def _fit_pls(gene_exp, imaging_data, n_components: int):
        return pls_regression(
            gene_exp,
            np.asarray(imaging_data, dtype=float).reshape(-1, 1),
            n_components=n_components,
            n_perm=0,
            n_boot=0,
        )

    @classmethod
    def set_coef(cls, data, gene_exp, var=None, n_components=None):
        max_components = min(15, int(np.asarray(data).shape[0]) - 1, int(np.asarray(gene_exp).shape[1]))
        if max_components < 1:
            raise ValueError("PLS requires at least two samples and one feature.")

        res = cls._fit_pls(gene_exp, data, n_components=max_components)
        explained_var = np.asarray(res.get("varexp"), dtype=float)
        cumulative_var = np.cumsum(explained_var)
        if n_components is None:
            cls.check_var(float(var))
            n_components = int(np.searchsorted(cumulative_var, float(var), side="left") + 1)
        elif n_components > max_components:
            raise ValueError(
                f"Requested n_components={n_components}, but only {max_components} components are available."
            )

        if var is None:
            var = float(cumulative_var[n_components - 1])
        return float(var), int(n_components), explained_var

    @property
    def p_val(self):
        return self._p_val

    @property
    def r2(self):
        return self._r2

    @p_val.setter
    def p_val(self, p_val):
        self._p_val = np.asarray(p_val, dtype=float)

    @r2.setter
    def r2(self, r2):
        self._r2 = np.asarray(r2, dtype=float)

    def boot_pls(self, imaging_data, permuted_imaging, gene_exp):  # pragma: no cover
        """Run permutation testing for the retained PLS components."""

        logger.info("Calculating PLS with permuted data")
        original = self._fit_pls(gene_exp, imaging_data, n_components=self.n_components)
        target_cumulative = np.cumsum(np.asarray(original.get("varexp"), dtype=float))[: self.n_components]
        self.r2[:] = target_cumulative

        permuted_cumulative = np.zeros((permuted_imaging.shape[1], self.n_components), dtype=float)
        for index in range(permuted_imaging.shape[1]):
            permuted = self._fit_pls(gene_exp, permuted_imaging[:, index], n_components=self.n_components)
            permuted_cumulative[index, :] = np.cumsum(np.asarray(permuted.get("varexp"), dtype=float))[
                : self.n_components
            ]
        self.p_val[:] = (np.sum(permuted_cumulative >= self.r2.reshape(1, -1), axis=0) + 1) / (
            permuted_imaging.shape[1] + 1
        )
        return
