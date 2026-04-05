from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

from ._logging import get_logger
from .genes import GeneResults
from .pls_backend import (
    PreparedPLS1,
    fit_prepared_pls1,
    fit_prepared_pls1_chunk,
    pls_regression,
    prepare_pls1,
)

logger = get_logger(__name__)


def _progress_step(total: int) -> int:
    """Choose a coarse logging interval for long permutation loops."""

    return max(1, total // 10)


def _log_permutation_progress(completed: int, total: int, *, step: int) -> None:
    """Emit periodic progress updates for long permutation loops."""

    if total < 1:
        return
    if completed == total or completed % step == 0:
        logger.info("Processed %d/%d PLS permutations.", completed, total)


def _chunk_bounds(total: int, n_jobs: int) -> list[tuple[int, int]]:
    """Split a permutation loop into a small number of larger work chunks."""

    n_chunks = min(total, max(1, n_jobs) * 4)
    chunk_size = max(1, (total + n_chunks - 1) // n_chunks)
    return [
        (start, min(total, start + chunk_size))
        for start in range(0, total, chunk_size)
    ]


class PLSAnalysis:
    """Perform PLS regression for the imaging transcriptomics workflows."""

    def __init__(
        self,
        imaging_data,
        gene_exp,
        n_components: int,
        var: float,
        n_iter: int = 1000,
        n_jobs: int = 1,
        *,
        store_gene_nulls: bool = True,
    ):
        """Fit the initial PLS model and prepare storage for permutation results."""

        self._prepared_gene_exp = prepare_pls1(gene_exp)
        self.var, self.n_components, self.components_var, self._initial_fit = self.set_coef(
            imaging_data,
            self._prepared_gene_exp,
            n_components=n_components,
            var=var,
        )
        self.n_iter = n_iter
        self.n_jobs = max(1, int(n_jobs))
        self._p_val = np.zeros(self.n_components, dtype=float)
        self._r2 = np.zeros(self.n_components, dtype=float)
        self.gene_results = GeneResults(
            "pls",
            n_components=self.n_components,
            n_iter=n_iter,
            n_genes=np.asarray(gene_exp).shape[1],
            store_weights=store_gene_nulls,
        )

    @staticmethod
    def check_var(var: float):
        """Validate a cumulative explained-variance target."""

        if var < 0:
            raise ValueError("The variance must be a positive number.")
        if var > 1:
            raise ValueError("The variance must be a number between 0 and 1.")
        return var

    @staticmethod
    def _fit_pls(
        gene_exp,
        imaging_data,
        n_components: int,
        *,
        return_x_scores: bool = True,
        return_x_weights: bool = True,
    ):
        """Fit one PLS-1 model using either a prepared or raw expression matrix."""

        if isinstance(gene_exp, PreparedPLS1):
            return fit_prepared_pls1(
                gene_exp,
                np.asarray(imaging_data, dtype=float).reshape(-1, 1),
                n_components=n_components,
                return_full=False,
                return_x_scores=return_x_scores,
                return_x_weights=return_x_weights,
            )
        return pls_regression(
            gene_exp,
            np.asarray(imaging_data, dtype=float).reshape(-1, 1),
            n_components=n_components,
            n_perm=0,
            n_boot=0,
        )

    @classmethod
    def set_coef(cls, data, gene_exp, var=None, n_components=None):
        """Choose the retained component count and initial explained variance."""

        n_features = int(gene_exp.X.shape[1]) if isinstance(gene_exp, PreparedPLS1) else int(np.asarray(gene_exp).shape[1])
        max_components = min(15, int(np.asarray(data).shape[0]) - 1, n_features)
        if max_components < 1:
            raise ValueError("PLS requires at least two samples and one feature.")

        res = cls._fit_pls(gene_exp, data, n_components=max_components, return_x_scores=True, return_x_weights=True)
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
        return float(var), int(n_components), explained_var, res

    @staticmethod
    def _truncate_fit(fit_result, n_components: int):
        """Keep only the fitted fields required by downstream code paths."""

        return {
            "x_scores": np.asarray(fit_result.get("x_scores"), dtype=float)[:, :n_components],
            "x_weights": np.asarray(fit_result.get("x_weights"), dtype=float)[:, :n_components],
            "varexp": np.asarray(fit_result.get("varexp"), dtype=float)[:n_components],
        }

    @property
    def p_val(self):
        """Nominal component p-values from permutation testing."""

        return self._p_val

    @property
    def r2(self):
        """Cumulative explained variance for the retained components."""

        return self._r2

    @p_val.setter
    def p_val(self, p_val):
        """Store component p-values as a float array."""

        self._p_val = np.asarray(p_val, dtype=float)

    @r2.setter
    def r2(self, r2):
        """Store cumulative explained variance as a float array."""

        self._r2 = np.asarray(r2, dtype=float)

    def boot_pls(
        self,
        imaging_data,
        permuted_imaging,
        gene_exp,
        *,
        scan_data=None,
        gene_labels=None,
    ):  # pragma: no cover
        """Run permutation testing for the retained PLS components.

        When ``scan_data`` and ``gene_labels`` are provided, this same loop also
        fills the gene-weight bootstrap used for downstream PLS gene tables.
        """

        total_permutations = int(permuted_imaging.shape[1])
        progress_step = _progress_step(total_permutations)
        logger.info(
            "Calculating PLS with permuted data (%d permutations, %d job%s).",
            total_permutations,
            self.n_jobs,
            "" if self.n_jobs == 1 else "s",
        )
        original = self._truncate_fit(self._initial_fit, self.n_components)
        target_cumulative = np.cumsum(np.asarray(original.get("varexp"), dtype=float))[: self.n_components]
        self.r2[:] = target_cumulative
        if scan_data is not None and gene_labels is not None:
            self.gene_results.results.prepare_from_fit(original, scan_data, gene_labels)

        permuted_cumulative = np.zeros((total_permutations, self.n_components), dtype=float)
        need_weights = scan_data is not None and gene_labels is not None
        if self.n_jobs == 1 or total_permutations <= 1:
            for index in range(total_permutations):
                permuted = self._fit_pls(
                    gene_exp,
                    permuted_imaging[:, index],
                    n_components=self.n_components,
                    return_x_scores=False,
                    return_x_weights=need_weights,
                )
                permuted_cumulative[index, :] = np.cumsum(np.asarray(permuted.get("varexp"), dtype=float))[
                    : self.n_components
                ]
                if need_weights:
                    self.gene_results.results.store_permuted_weights(index, permuted.get("x_weights"))
                _log_permutation_progress(index + 1, total_permutations, step=progress_step)
        else:
            max_workers = min(self.n_jobs, total_permutations)

            def _fit_chunk(bounds: tuple[int, int]):
                start, end = bounds
                if isinstance(gene_exp, PreparedPLS1):
                    block_fit = fit_prepared_pls1_chunk(
                        gene_exp,
                        permuted_imaging[:, start:end],
                        n_components=self.n_components,
                        return_x_weights=need_weights,
                    )
                    cumulative_block = np.cumsum(np.asarray(block_fit["varexp"], dtype=float), axis=1)[
                        :, : self.n_components
                    ]
                    weight_block = (
                        np.asarray(block_fit["x_weights"], dtype=float)
                        if need_weights and "x_weights" in block_fit
                        else None
                    )
                    return start, end, cumulative_block, weight_block

                cumulative_block = np.zeros((end - start, self.n_components), dtype=float)
                weight_block = [] if need_weights else None
                for local_index, index in enumerate(range(start, end)):
                    permuted = self._fit_pls(
                        gene_exp,
                        permuted_imaging[:, index],
                        n_components=self.n_components,
                        return_x_scores=False,
                        return_x_weights=need_weights,
                    )
                    cumulative_block[local_index, :] = np.cumsum(
                        np.asarray(permuted.get("varexp"), dtype=float)
                    )[: self.n_components]
                    if need_weights:
                        weight_block.append(np.asarray(permuted.get("x_weights"), dtype=float))
                return start, end, cumulative_block, weight_block

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [
                    executor.submit(_fit_chunk, bounds)
                    for bounds in _chunk_bounds(total_permutations, max_workers)
                ]
                completed = 0
                for future in as_completed(futures):
                    start, end, cumulative_block, weight_block = future.result()
                    permuted_cumulative[start:end, :] = cumulative_block
                    if need_weights and weight_block is not None:
                        for offset, x_weights in enumerate(weight_block):
                            self.gene_results.results.store_permuted_weights(start + offset, x_weights)
                    completed += end - start
                    _log_permutation_progress(completed, total_permutations, step=progress_step)
        logger.info("Finished PLS permutation fits.")
        self.p_val[:] = (np.sum(permuted_cumulative >= self.r2.reshape(1, -1), axis=0) + 1) / (
            total_permutations + 1
        )
        return
