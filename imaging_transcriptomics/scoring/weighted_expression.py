"""Numerical helpers for PTRS-style weighted regional gene-expression scores."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..exceptions import ConfigurationError
from ..models import GEDARDirection, WeightNormalization


@dataclass(frozen=True)
class WeightedProjection:
    """The weights used for one projection plus the resulting regional scores."""

    weights: np.ndarray
    selected: np.ndarray
    score: np.ndarray
    score_z: np.ndarray


def direction_modes(direction: GEDARDirection) -> tuple[str, ...]:
    """Return the concrete direction modes implied by one GEDAR option."""

    if direction == "split":
        return ("up", "down")
    if direction in {"combined", "up", "down"}:
        return (direction,)
    raise ConfigurationError("direction must be one of 'combined', 'up', 'down', or 'split'.")


def direction_mask(weights: np.ndarray, direction: str) -> np.ndarray:
    """Build the initial boolean mask for the requested weight direction."""

    if direction == "combined":
        return np.ones(weights.shape[0], dtype=bool)
    if direction == "up":
        return weights > 0
    if direction == "down":
        return weights < 0
    raise ConfigurationError("direction must be one of 'combined', 'up', or 'down'.")


def normalize_weight_vector(weights: np.ndarray, mode: WeightNormalization) -> np.ndarray:
    """Normalize one selected weight vector with the configured mode."""

    vector = np.asarray(weights, dtype=float).reshape(-1)
    if mode == "none":
        return vector
    if mode == "zscore":
        centered = vector - vector.mean()
        scale = vector.std(ddof=1)
        if not np.isfinite(scale) or scale == 0:
            return centered
        return centered / scale
    if mode == "unit":
        norm = float(np.linalg.norm(vector))
        if not np.isfinite(norm) or norm == 0:
            return vector
        return vector / norm
    raise ConfigurationError("normalize_weights must be one of 'none', 'zscore', or 'unit'.")


def standardize_vector(values: np.ndarray) -> np.ndarray:
    """Return a z-scored copy of one regional score vector."""

    vector = np.asarray(values, dtype=float).reshape(-1)
    centered = vector - vector.mean()
    scale = vector.std(ddof=1)
    if not np.isfinite(scale) or scale == 0:
        return centered
    return centered / scale


def projection_weights(weights: np.ndarray, direction: str) -> np.ndarray:
    """Return the effective weights used by the PTRS-style GEDAR average."""

    vector = np.asarray(weights, dtype=float).reshape(-1)
    if direction == "combined":
        return vector
    if direction in {"up", "down"}:
        return np.abs(vector)
    raise ConfigurationError("direction must be one of 'combined', 'up', or 'down'.")


def score_regions(expression: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Compute a PTRS-style weighted regional average from expression and weights."""

    weight_vector = np.asarray(weights, dtype=float).reshape(-1)
    denominator = float(weight_vector.sum())
    if not np.isfinite(denominator) or np.isclose(denominator, 0.0):
        raise ConfigurationError(
            "The selected GEDAR weights sum to zero, so the weighted average is undefined. "
            "Use normalize_weights='none' or 'unit', or choose a one-sided direction."
        )
    return (np.asarray(expression, dtype=float) @ weight_vector) / denominator


def project_weighted_expression(
    expression_matrix: np.ndarray,
    input_weights: np.ndarray,
    selected_mask: np.ndarray,
    *,
    direction: str,
    normalize_weights: WeightNormalization,
    n_regions: int,
) -> WeightedProjection:
    """Project selected gene weights onto regional expression values."""

    selected = np.asarray(selected_mask, dtype=bool).reshape(-1)
    weights = np.zeros(selected.shape[0], dtype=float)
    if not np.any(selected):
        nan_vector = np.full(n_regions, np.nan, dtype=float)
        return WeightedProjection(weights=weights, selected=selected, score=nan_vector, score_z=nan_vector.copy())

    normalized = normalize_weight_vector(
        np.asarray(input_weights, dtype=float)[selected],
        normalize_weights,
    )
    projected = projection_weights(normalized, direction)
    weights[selected] = projected
    score = score_regions(np.asarray(expression_matrix, dtype=float)[:, selected], projected)
    return WeightedProjection(
        weights=weights,
        selected=selected,
        score=score,
        score_z=standardize_vector(score),
    )
