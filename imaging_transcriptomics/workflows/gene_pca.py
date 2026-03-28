"""Workflow implementation for atlas-based gene-expression PCA."""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from ..exceptions import ConfigurationError
from ..gene_expression import load_brain_gene_symbols
from ..models import AtlasSelection, GenePCAResult, HemisphereMode, RegionScope
from ..preprocessing.gene_lists import (
    load_gene_list,
    partition_brain_filtered_genes,
    resolve_selected_genes,
)


def _standardize_columns(matrix: np.ndarray) -> np.ndarray:
    """Z-score each gene across regions before PCA."""

    numeric = np.asarray(matrix, dtype=float)
    finite_mask = np.isfinite(numeric)
    safe = np.where(finite_mask, numeric, 0.0)
    counts = finite_mask.sum(axis=0, keepdims=True)
    means = np.divide(
        safe.sum(axis=0, keepdims=True),
        counts,
        out=np.zeros((1, numeric.shape[1]), dtype=float),
        where=counts > 0,
    )
    centered = np.where(finite_mask, numeric - means, 0.0)
    centered_ss = np.sum(centered * centered, axis=0, keepdims=True)
    scale = np.sqrt(
        np.divide(
            centered_ss,
            np.maximum(counts - 1, 1),
            out=np.zeros((1, numeric.shape[1]), dtype=float),
            where=counts > 1,
        )
    )
    scale[scale == 0] = 1.0
    standardized = centered / scale
    standardized[~np.isfinite(standardized)] = 0.0
    return standardized


def _fit_pca(matrix: np.ndarray, n_components: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fit PCA with SVD and return scores, loadings, and explained variance ratio."""

    standardized = _standardize_columns(np.asarray(matrix, dtype=float))
    try:
        u, singular_values, vt = np.linalg.svd(standardized, full_matrices=False)
    except np.linalg.LinAlgError as exc:
        raise ConfigurationError(
            "PCA failed to converge on the selected atlas expression matrix. "
            "This usually means the selected genes contain too many invalid values. "
            "Try a different gene list or inspect the atlas expression inputs."
        ) from exc

    scores = u[:, :n_components] * singular_values[:n_components]
    loadings = vt[:n_components, :].T
    if standardized.shape[0] > 1:
        explained_variance = (singular_values[:n_components] ** 2) / (standardized.shape[0] - 1)
        total_variance = (singular_values**2).sum() / (standardized.shape[0] - 1)
    else:
        explained_variance = singular_values[:n_components] ** 2
        total_variance = (singular_values**2).sum()
    explained_ratio = np.zeros(n_components, dtype=float) if total_variance == 0 else explained_variance / total_variance
    return scores, loadings, explained_ratio


def _component_columns(prefix: str, n_components: int) -> list[str]:
    """Create predictable component column names such as PC1 and PC2."""

    return [f"{prefix}{index}" for index in range(1, n_components + 1)]


def _regional_scores_frame(labels: pd.DataFrame, scores: np.ndarray) -> pd.DataFrame:
    """Build the regional score table for all retained principal components."""

    return labels.reset_index(drop=True).assign(
        **dict(zip(_component_columns("PC", scores.shape[1]), scores.T, strict=False))
    )


def _gene_loadings_frame(matched_genes: list[str], loadings: np.ndarray) -> pd.DataFrame:
    """Build the per-gene loading table for all retained components."""

    return pd.DataFrame(
        {
            "gene": matched_genes,
            **dict(zip(_component_columns("PC", loadings.shape[1]), loadings.T, strict=False)),
        }
    )


def _variance_frame(explained_ratio: np.ndarray) -> pd.DataFrame:
    """Build the explained-variance summary table for the PCA run."""

    cumulative = np.cumsum(explained_ratio)
    return pd.DataFrame(
        {
            "component": np.arange(1, explained_ratio.shape[0] + 1),
            "variance_explained": explained_ratio,
            "cumulative_variance": cumulative,
        }
    )


def run_gene_pca(
    genes: str | Path | Iterable[str],
    *,
    atlas: str = "dk",
    hemisphere: HemisphereMode = "left",
    regions: RegionScope = "all",
    n_components: int = 3,
    output_dir: str | Path | None = None,
    select_atlas_data_fn: Callable[..., AtlasSelection],
    write_result_bundle_fn: Callable[[GenePCAResult, Path], None] | None = None,
) -> GenePCAResult:
    """Run PCA on atlas expression after filtering to a selected gene list."""

    requested_genes = load_gene_list(genes)
    if not requested_genes:
        raise ConfigurationError("No gene symbols were provided.")
    if int(n_components) < 1:
        raise ConfigurationError("n_components must be at least 1.")

    filtered_request, brain_filtered_genes = partition_brain_filtered_genes(
        requested_genes,
        load_brain_gene_symbols(),
    )
    if not filtered_request:
        raise ConfigurationError("No genes remained after applying the packaged AHPA brain-gene filter.")

    selection = select_atlas_data_fn(atlas=atlas, hemisphere=hemisphere, regions=regions)
    available_genes = selection.expression.columns[2:].to_numpy(dtype=str, copy=False)
    matched_genes, missing_genes = resolve_selected_genes(filtered_request, available_genes)
    if not matched_genes:
        raise ConfigurationError(
            f"None of the requested genes were found in atlas '{selection.atlas.id}'."
        )

    max_components = min(int(n_components), selection.n_regions, len(matched_genes))
    gene_matrix = selection.expression.loc[:, matched_genes].to_numpy(dtype=float, copy=True)
    scores, loadings, explained_ratio = _fit_pca(gene_matrix, max_components)
    result = GenePCAResult(
        atlas_id=selection.atlas.id,
        atlas_label=selection.atlas.label,
        hemisphere=selection.hemisphere,
        regions=selection.regions,
        requested_genes=tuple(requested_genes),
        regional_scores=_regional_scores_frame(selection.labels, scores),
        gene_loadings=_gene_loadings_frame(matched_genes, loadings),
        variance_table=_variance_frame(explained_ratio),
        matched_genes=tuple(matched_genes),
        brain_filtered_genes=tuple(brain_filtered_genes),
        missing_genes=tuple(missing_genes),
        output_dir=None if output_dir is None else Path(output_dir),
    )
    if output_dir is not None and write_result_bundle_fn is not None:
        write_result_bundle_fn(result, Path(output_dir))
    return result
