from __future__ import annotations

import logging
from contextlib import contextmanager

import numpy as np
import pandas as pd


def make_prerank_table(genes, scores) -> pd.DataFrame:
    """Return a deterministic prerank table with duplicate scores broken by tiny offsets.

    GSEA libraries warn when tied ranking scores appear because the internal
    sort order for those genes can become arbitrary. We preserve the incoming
    gene order and add epsilon-scale offsets only within exact-tie groups.
    """

    gene_array = np.asarray(genes, dtype=object).reshape(-1)
    score_array = np.asarray(scores, dtype=float).reshape(-1)
    if gene_array.shape[0] != score_array.shape[0]:
        raise ValueError("Genes and scores must have the same length.")

    adjusted = score_array.copy()
    unique_vals, inverse, counts = np.unique(score_array, return_inverse=True, return_counts=True)
    for group_index, count in enumerate(counts):
        if count <= 1:
            continue
        member_idx = np.flatnonzero(inverse == group_index)
        base = float(unique_vals[group_index])
        left_gap = np.inf if group_index == 0 else base - float(unique_vals[group_index - 1])
        right_gap = np.inf if group_index == len(unique_vals) - 1 else float(unique_vals[group_index + 1]) - base
        min_gap = min(left_gap, right_gap)
        eps_step = np.finfo(float).eps * max(1.0, abs(base))
        if np.isfinite(min_gap) and min_gap > 0:
            step = min(eps_step, (min_gap / (count + 1)) * 0.5)
        else:
            step = eps_step
        offsets = (np.arange(count, dtype=float) - ((count - 1) / 2.0)) * step
        adjusted[member_idx] = base + offsets

    return pd.DataFrame({"gene": gene_array, "score": adjusted})


class _SuppressDuplicatePrerankWarnings(logging.Filter):
    """Drop the noisy GSEApy warning about duplicated preranked scores."""

    def filter(self, record: logging.LogRecord) -> bool:
        return "Duplicated values found in preranked stats" not in record.getMessage()


@contextmanager
def _suppress_gseapy_duplicate_prerank_warnings():
    """Temporarily patch GSEApy logger setup to filter duplicate-score warnings."""

    try:
        import gseapy.base as gseapy_base
        import gseapy.utils as gseapy_utils
    except ImportError:
        yield
        return

    original_base_log_init = gseapy_base.log_init
    original_utils_log_init = gseapy_utils.log_init

    def quiet_log_init(name, log_level=logging.INFO, filename=None):
        logger = original_base_log_init(name, log_level=log_level, filename=filename)
        if not any(isinstance(item, _SuppressDuplicatePrerankWarnings) for item in logger.filters):
            logger.addFilter(_SuppressDuplicatePrerankWarnings())
        return logger

    gseapy_base.log_init = quiet_log_init
    gseapy_utils.log_init = quiet_log_init
    try:
        yield
    finally:
        gseapy_base.log_init = original_base_log_init
        gseapy_utils.log_init = original_utils_log_init


def run_prerank(gseapy, rnk, gene_set, **kwargs):
    """Run GSEApy prerank while suppressing duplicate-score warning noise."""

    with _suppress_gseapy_duplicate_prerank_warnings():
        return gseapy.prerank(rnk, gene_set, **kwargs)


def _validate_es_inputs(es, esnull) -> tuple[np.ndarray, np.ndarray]:
    """Validate observed and null enrichment score arrays and coerce shapes."""

    es_array = np.asarray(es, dtype=float).reshape(-1)
    esnull_array = np.asarray(esnull, dtype=float)
    if esnull_array.ndim != 2:
        raise ValueError("Enrichment nulls must be a 2D array.")
    if esnull_array.shape[0] != es_array.shape[0]:
        raise ValueError("Observed ES and null ES must agree on the number of terms.")
    return es_array, esnull_array


def _same_sign_null_means(es_array: np.ndarray, esnull_array: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return per-term positive and negative null ES means."""

    pos_mask = esnull_array >= 0
    neg_mask = esnull_array < 0

    pos_mean = np.divide(
        np.where(pos_mask, esnull_array, 0.0).sum(axis=1),
        pos_mask.sum(axis=1),
        out=np.full(es_array.shape, np.nan, dtype=float),
        where=pos_mask.any(axis=1),
    )
    neg_mean = np.divide(
        np.where(neg_mask, esnull_array, 0.0).sum(axis=1),
        neg_mask.sum(axis=1),
        out=np.full(es_array.shape, np.nan, dtype=float),
        where=neg_mask.any(axis=1),
    )
    return pos_mean, neg_mean


def normalize_enrichment_scores(es, esnull) -> np.ndarray:
    """Normalize observed ES values using same-sign null means, matching GSEApy."""

    es_array, esnull_array = _validate_es_inputs(es, esnull)
    pos_mean, neg_mean = _same_sign_null_means(es_array, esnull_array)

    nes = np.full(es_array.shape, np.nan, dtype=float)
    pos_terms = es_array >= 0
    neg_terms = ~pos_terms

    if np.any(pos_terms):
        nes[pos_terms] = np.divide(
            es_array[pos_terms],
            pos_mean[pos_terms],
            out=np.full(pos_terms.sum(), np.nan, dtype=float),
            where=np.isfinite(pos_mean[pos_terms]) & (pos_mean[pos_terms] != 0),
        )
    if np.any(neg_terms):
        nes[neg_terms] = np.divide(
            -es_array[neg_terms],
            neg_mean[neg_terms],
            out=np.full(neg_terms.sum(), np.nan, dtype=float),
            where=np.isfinite(neg_mean[neg_terms]) & (neg_mean[neg_terms] != 0),
        )

    return nes


def normalize_enrichment_nulls(es, esnull) -> np.ndarray:
    """Normalize null ES values using same-sign means, matching GSEApy logic."""

    es_array, esnull_array = _validate_es_inputs(es, esnull)
    pos_mean, neg_mean = _same_sign_null_means(es_array, esnull_array)

    return np.where(
        esnull_array >= 0,
        np.divide(
            esnull_array,
            pos_mean[:, np.newaxis],
            out=np.full(esnull_array.shape, np.nan, dtype=float),
            where=np.isfinite(pos_mean[:, np.newaxis]) & (pos_mean[:, np.newaxis] != 0),
        ),
        np.divide(
            -esnull_array,
            neg_mean[:, np.newaxis],
            out=np.full(esnull_array.shape, np.nan, dtype=float),
            where=np.isfinite(neg_mean[:, np.newaxis]) & (neg_mean[:, np.newaxis] != 0),
        ),
    )


def nominal_pvalues_from_nulls(es, esnull) -> np.ndarray:
    """Return sign-aware nominal p-values from external enrichment nulls."""

    es_array, esnull_array = _validate_es_inputs(es, esnull)
    pvals = np.zeros(es_array.shape, dtype=float)
    n_perm = esnull_array.shape[1]
    for index, observed in enumerate(es_array):
        count = (
            np.sum(esnull_array[index, :] >= observed)
            if observed >= 0
            else np.sum(esnull_array[index, :] <= observed)
        )
        pvals[index] = (count + 1) / (n_perm + 1)
    return pvals


def gsea_style_fdr(nes, nesnull) -> np.ndarray:
    """Return GSEA-style FDR q-values from observed and null NES arrays."""

    nes_array = np.asarray(nes, dtype=float).reshape(-1)
    nesnull_array = np.asarray(nesnull, dtype=float)
    if nesnull_array.ndim != 2:
        raise ValueError("Normalized null enrichment scores must be a 2D array.")
    if nesnull_array.shape[0] != nes_array.shape[0]:
        raise ValueError("Observed NES and null NES must agree on the number of terms.")

    finite_nes = nes_array[np.isfinite(nes_array)]
    finite_nesnull = nesnull_array[np.isfinite(nesnull_array)]
    if finite_nes.size == 0 or finite_nesnull.size == 0:
        return np.full(nes_array.shape, np.nan, dtype=float)

    sorted_null = np.sort(finite_nesnull)
    sorted_obs = np.sort(finite_nes)
    fdr = np.full(nes_array.shape, np.nan, dtype=float)

    for index, value in enumerate(nes_array):
        if not np.isfinite(value):
            continue
        if value >= 0:
            all_pos = int(sorted_null.size - np.searchsorted(sorted_null, 0, side="left"))
            all_higher_and_pos = int(sorted_null.size - np.searchsorted(sorted_null, value, side="left"))
            obs_pos = int(sorted_obs.size - np.searchsorted(sorted_obs, 0, side="left"))
            obs_higher_and_pos = int(sorted_obs.size - np.searchsorted(sorted_obs, value, side="left"))
        else:
            all_pos = int(np.searchsorted(sorted_null, 0, side="left"))
            all_higher_and_pos = int(np.searchsorted(sorted_null, value, side="right"))
            obs_pos = int(np.searchsorted(sorted_obs, 0, side="left"))
            obs_higher_and_pos = int(np.searchsorted(sorted_obs, value, side="right"))

        if all_pos == 0 or obs_pos == 0 or obs_higher_and_pos == 0:
            fdr[index] = np.nan
            continue

        pi_norm = all_higher_and_pos / float(all_pos)
        pi_obs = obs_higher_and_pos / float(obs_pos)
        value_fdr = pi_norm / pi_obs if pi_obs != 0 else np.nan
        fdr[index] = min(value_fdr, 1.0) if np.isfinite(value_fdr) else np.nan

    return fdr
