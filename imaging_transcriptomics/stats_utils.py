from __future__ import annotations

import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """Return Benjamini-Hochberg adjusted p-values."""

    return multipletests(np.asarray(p_values, dtype=float), method="fdr_bh", is_sorted=False)[1]


def two_sided_z_pvalues(z_scores: np.ndarray) -> np.ndarray:
    """Return two-sided z-test p-values."""

    return 2 * norm.sf(np.abs(np.asarray(z_scores, dtype=float)))


def empirical_signed_pvalues(observed: np.ndarray, null_values: np.ndarray) -> np.ndarray:
    """Return nominal p-values using sign-aware permutation comparisons."""

    obs = np.asarray(observed, dtype=float).reshape(-1)
    nulls = np.asarray(null_values, dtype=float)
    if nulls.ndim != 2 or nulls.shape[0] != obs.shape[0]:
        raise ValueError("Null values must be a 2D array with one row per observed statistic.")
    pos_counts = np.sum(nulls >= obs.reshape(-1, 1), axis=1)
    neg_counts = np.sum(nulls <= obs.reshape(-1, 1), axis=1)
    counts = np.where(obs >= 0, pos_counts, neg_counts)
    return (counts + 1) / (nulls.shape[1] + 1)


def max_t_fwer_abs(observed: np.ndarray, null_values: np.ndarray) -> np.ndarray:
    """Return maxT family-wise corrected p-values using absolute test statistics."""

    obs = np.abs(np.asarray(observed, dtype=float).reshape(-1))
    nulls = np.asarray(null_values, dtype=float)
    if nulls.ndim != 2:
        raise ValueError("Null values must be a 2D array.")
    perm_max_abs = np.max(np.abs(nulls), axis=0)
    counts = np.sum(perm_max_abs.reshape(1, -1) >= obs.reshape(-1, 1), axis=1)
    return (counts + 1) / (nulls.shape[1] + 1)


def minimum_bh_resolution(n_tests: int, n_permutations: int) -> float:
    """Return the smallest achievable BH-adjusted p-value for a permutation grid."""

    return float(n_tests) / float(n_permutations + 1)
