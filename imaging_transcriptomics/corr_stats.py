from __future__ import annotations

import numpy as np


def rank_vector(values: np.ndarray) -> np.ndarray:
    """Return stable zero-based ranks for a one-dimensional array."""

    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.shape[0], dtype=float)
    ranks[order] = np.arange(values.shape[0], dtype=float)
    return ranks


def _centered_rank_vector(values: np.ndarray) -> np.ndarray:
    """Return integer-centered stable ranks for a one-dimensional array."""

    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.shape[0], dtype=np.int64)
    ranks[order] = np.arange(values.shape[0], dtype=np.int64)
    return (ranks * 2) - (values.shape[0] - 1)


def rank_columns(matrix: np.ndarray) -> np.ndarray:
    """Return stable zero-based ranks for each column of a matrix."""

    order = np.argsort(matrix, axis=0, kind="mergesort")
    ranks = np.empty(order.shape, dtype=float)
    template = np.arange(matrix.shape[0], dtype=float)[:, None]
    np.put_along_axis(ranks, order, template, axis=0)
    return ranks


def _centered_rank_columns(matrix: np.ndarray) -> np.ndarray:
    """Return integer-centered stable ranks for each column of a matrix."""

    order = np.argsort(matrix, axis=0, kind="mergesort")
    ranks = np.empty(order.shape, dtype=np.int64)
    template = np.arange(matrix.shape[0], dtype=np.int64)[:, None]
    np.put_along_axis(ranks, order, template, axis=0)
    return (ranks * 2) - (matrix.shape[0] - 1)


def standardize_columns(matrix: np.ndarray) -> np.ndarray:
    """Center and scale each column with sample standard deviation."""

    centered = matrix - matrix.mean(axis=0, keepdims=True)
    scale = matrix.std(axis=0, ddof=1, keepdims=True)
    scale[scale == 0] = 1.0
    return centered / scale


def standardize_vector(values: np.ndarray) -> np.ndarray:
    """Center and scale a vector with sample standard deviation."""

    centered = values - values.mean()
    scale = values.std(ddof=1)
    if scale == 0:
        return centered
    return centered / scale


def ranked_gene_expression(gene_exp: np.ndarray) -> np.ndarray:
    """Pre-rank atlas gene expression for deterministic Spearman correlation."""

    return _centered_rank_columns(np.asarray(gene_exp, dtype=float))


def spearman_correlation_matrix(imaging_data: np.ndarray, ranked_genes: np.ndarray) -> np.ndarray:
    """Compute Spearman correlation between one imaging vector and all genes."""

    ranked_img = _centered_rank_vector(np.asarray(imaging_data, dtype=float).reshape(-1))
    n = ranked_img.shape[0]
    denom = max(n * (n * n - 1), 1)
    return (3.0 * (np.asarray(ranked_genes, dtype=np.int64).T @ ranked_img.reshape(-1, 1))) / denom


def spearman_correlation_bootstrap(permuted_imaging: np.ndarray, ranked_genes: np.ndarray) -> np.ndarray:
    """Compute Spearman correlations for all permuted imaging vectors at once."""

    ranked_perm = _centered_rank_columns(np.asarray(permuted_imaging, dtype=float))
    n = ranked_perm.shape[0]
    denom = max(n * (n * n - 1), 1)
    return (3.0 * (np.asarray(ranked_genes, dtype=np.int64).T @ ranked_perm)) / denom
