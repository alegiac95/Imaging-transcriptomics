from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .exceptions import ConfigurationError
from .gene_expression import select_atlas_data
from .models import GenePCAResult
from .serialization import write_result_bundle


def _dedupe_preserve_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        item = str(value).strip()
        if not item:
            continue
        key = item.upper()
        if key in seen:
            continue
        seen.add(key)
        out.append(item)
    return out


def _parse_gene_tokens(text: str) -> list[str]:
    tokens: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        normalized = line.replace(",", " ").replace("\t", " ")
        tokens.extend(part for part in normalized.split() if part)
    return _dedupe_preserve_order(tokens)


def load_gene_list(genes: str | Path | Iterable[str]) -> list[str]:
    if isinstance(genes, (list, tuple, set, np.ndarray, pd.Series)):
        return _dedupe_preserve_order(str(item) for item in genes)
    path = Path(str(genes)).expanduser()
    if path.exists() and path.is_file():
        return _parse_gene_tokens(path.read_text())
    return _parse_gene_tokens(str(genes))


def _standardize_columns(matrix: np.ndarray) -> np.ndarray:
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    scale = matrix.std(axis=0, ddof=1, keepdims=True)
    scale[~np.isfinite(scale) | (scale == 0)] = 1.0
    return centered / scale


def _fit_pca(matrix: np.ndarray, n_components: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    standardized = _standardize_columns(np.asarray(matrix, dtype=float))
    u, singular_values, vt = np.linalg.svd(standardized, full_matrices=False)
    scores = u[:, :n_components] * singular_values[:n_components]
    loadings = vt[:n_components, :].T
    if standardized.shape[0] > 1:
        explained_variance = (singular_values[:n_components] ** 2) / (standardized.shape[0] - 1)
        total_variance = (singular_values ** 2).sum() / (standardized.shape[0] - 1)
    else:
        explained_variance = singular_values[:n_components] ** 2
        total_variance = (singular_values ** 2).sum()
    if total_variance == 0:
        explained_ratio = np.zeros(n_components, dtype=float)
    else:
        explained_ratio = explained_variance / total_variance
    return scores, loadings, explained_ratio


def _resolve_selected_genes(requested_genes: list[str], available_genes: np.ndarray) -> tuple[list[str], list[str]]:
    lookup: dict[str, str] = {}
    for gene in available_genes.astype(str).tolist():
        lookup.setdefault(gene.upper(), gene)
    matched: list[str] = []
    missing: list[str] = []
    for gene in requested_genes:
        canonical = lookup.get(gene.upper())
        if canonical is None:
            missing.append(gene)
        else:
            matched.append(canonical)
    return matched, missing


def _component_columns(prefix: str, n_components: int) -> list[str]:
    return [f"{prefix}{index}" for index in range(1, n_components + 1)]


def _regional_scores_frame(labels: pd.DataFrame, scores: np.ndarray) -> pd.DataFrame:
    return labels.reset_index(drop=True).assign(**dict(zip(_component_columns("PC", scores.shape[1]), scores.T, strict=False)))


def _gene_loadings_frame(matched_genes: list[str], loadings: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene": matched_genes,
            **dict(zip(_component_columns("PC", loadings.shape[1]), loadings.T, strict=False)),
        }
    )


def _variance_frame(explained_ratio: np.ndarray) -> pd.DataFrame:
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
    hemisphere: str = "left",
    regions: str = "all",
    n_components: int = 3,
    output_dir: str | Path | None = None,
) -> GenePCAResult:
    """Run PCA on atlas expression after filtering to a selected gene list."""

    requested_genes = load_gene_list(genes)
    if not requested_genes:
        raise ConfigurationError("No gene symbols were provided.")
    if int(n_components) < 1:
        raise ConfigurationError("n_components must be at least 1.")

    selection = select_atlas_data(atlas=atlas, hemisphere=hemisphere, regions=regions)
    available_genes = selection.expression.columns[2:].to_numpy(dtype=str, copy=False)
    matched_genes, missing_genes = _resolve_selected_genes(requested_genes, available_genes)
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
        missing_genes=tuple(missing_genes),
        output_dir=None if output_dir is None else Path(output_dir),
    )
    if output_dir is not None:
        write_result_bundle(result, Path(output_dir))
    return result
