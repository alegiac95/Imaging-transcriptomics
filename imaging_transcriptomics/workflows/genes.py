"""Workflow implementations for single-gene atlas queries."""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from ..corr_stats import ranked_gene_expression, spearman_correlation_matrix
from ..exceptions import ConfigurationError
from ..models import (
    AtlasSelection,
    GeneQueryResult,
    HemisphereMode,
    RegionScope,
)
from ..stats_utils import bh_fdr


def _resolve_gene(selection: AtlasSelection, gene: str) -> str:
    """Resolve one user-provided gene symbol against atlas columns."""

    query = str(gene).strip().upper()
    if not query:
        raise ConfigurationError("A non-empty gene symbol is required.")
    atlas_genes = selection.expression.columns[2:].astype(str)
    gene_lookup = {symbol.upper(): symbol for symbol in atlas_genes}
    matched = gene_lookup.get(query)
    if matched is None:
        raise ConfigurationError(f"Gene '{query}' was not found in atlas '{selection.atlas.id}'.")
    return str(matched)


def _regional_expression_frame(selection: AtlasSelection, gene: str, *, zscore_expression: bool) -> pd.DataFrame:
    """Return one atlas-aligned regional expression table for a selected gene."""

    value_column = "expression_z" if zscore_expression else "expression"
    return selection.labels.reset_index(drop=True).assign(
        **{value_column: selection.expression[gene].to_numpy(dtype=float, copy=True)}
    )


def _coexpression_table(
    selection: AtlasSelection,
    gene: str,
    *,
    fdr_threshold: float,
    top_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return full, top-positive, and top-negative co-expression tables."""

    if selection.n_regions < 4:
        raise ConfigurationError("At least four atlas regions are required to estimate co-expression.")

    atlas_genes = selection.expression.columns[2:].astype(str).to_numpy(dtype=str, copy=True)
    target_idx = int(np.flatnonzero(atlas_genes == gene)[0])
    expression_matrix = selection.expression.iloc[:, 2:].to_numpy(dtype=float, copy=True)
    correlations = spearman_correlation_matrix(
        expression_matrix[:, target_idx],
        ranked_gene_expression(expression_matrix),
    ).reshape(-1)

    n_regions = selection.n_regions
    df = max(n_regions - 2, 1)
    clipped = np.clip(correlations, -0.999999, 0.999999)
    with np.errstate(divide="ignore", invalid="ignore"):
        t_stat = clipped * np.sqrt(df / np.maximum(1.0 - (clipped**2), np.finfo(float).eps))
    p_values = 2.0 * student_t.sf(np.abs(t_stat), df)
    p_values[~np.isfinite(p_values)] = 1.0

    table = pd.DataFrame(
        {
            "gene": atlas_genes,
            "score": correlations,
            "p": p_values,
        }
    )
    table = table.loc[table["gene"] != gene].copy()
    table["fdr"] = bh_fdr(table["p"].to_numpy(dtype=float))
    scores = table["score"].to_numpy(dtype=float)
    table["significant"] = table["fdr"].to_numpy(dtype=float) <= float(fdr_threshold)
    table["selected"] = table["significant"] & (scores > 0)
    table["selected_negative"] = table["significant"] & (scores < 0)
    table = table.sort_values(
        ["selected", "selected_negative", "score", "gene"],
        ascending=[False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    top_positive = table.loc[table["selected"]].head(int(top_n)).reset_index(drop=True)
    top_negative = (
        table.loc[table["selected_negative"]]
        .sort_values(["score", "gene"], ascending=[True, True], kind="mergesort")
        .head(int(top_n))
        .reset_index(drop=True)
    )
    return table, top_positive, top_negative


def _coexpression_matrix(
    selection: AtlasSelection,
    gene: str,
    top_positive: pd.DataFrame,
    top_negative: pd.DataFrame,
) -> pd.DataFrame:
    """Return a seed-plus-top-genes Spearman co-expression matrix."""

    ordered_genes = [
        gene,
        *top_positive["gene"].astype(str).tolist(),
        *top_negative["gene"].astype(str).tolist(),
    ]
    seen: set[str] = set()
    unique_genes = [symbol for symbol in ordered_genes if not (symbol in seen or seen.add(symbol))]
    matrix = selection.expression.loc[:, unique_genes].corr(method="spearman")
    return matrix.loc[unique_genes, unique_genes]


def run_gene_query(
    gene: str,
    *,
    atlas: str = "dk",
    hemisphere: HemisphereMode = "left",
    regions: RegionScope = "default",
    zscore_expression: bool = True,
    top_n: int = 25,
    fdr_threshold: float = 0.05,
    output_dir: str | Path | None = None,
    select_atlas_data_fn: Callable[..., AtlasSelection],
    write_result_bundle_fn: Callable[[GeneQueryResult, Path], None] | None = None,
) -> GeneQueryResult:
    """Return one gene's regional expression together with its top co-expressed genes."""

    if int(top_n) < 1:
        raise ConfigurationError("top_n must be at least 1.")
    if not (0.0 < float(fdr_threshold) <= 1.0):
        raise ConfigurationError("fdr_threshold must be between 0 and 1.")

    selection = select_atlas_data_fn(
        atlas=atlas,
        hemisphere=hemisphere,
        regions=regions,
        zscore_expression=zscore_expression,
    )
    matched_gene = _resolve_gene(selection, gene)
    gene_table, top_positive, top_negative = _coexpression_table(
        selection,
        matched_gene,
        fdr_threshold=float(fdr_threshold),
        top_n=int(top_n),
    )
    coexpression_matrix = _coexpression_matrix(selection, matched_gene, top_positive, top_negative)
    result = GeneQueryResult(
        atlas_id=selection.atlas.id,
        atlas_label=selection.atlas.label,
        hemisphere=selection.hemisphere,
        regions=selection.regions,
        gene=matched_gene,
        zscore_expression=bool(zscore_expression),
        regional_values=_regional_expression_frame(selection, matched_gene, zscore_expression=bool(zscore_expression)),
        gene_table=gene_table,
        coexpressed_genes=top_positive,
        anticorrelated_genes=top_negative,
        coexpression_matrix=coexpression_matrix,
        top_n=int(top_n),
        fdr_threshold=float(fdr_threshold),
        output_dir=None if output_dir is None else Path(output_dir),
    )
    if output_dir is not None and write_result_bundle_fn is not None:
        write_result_bundle_fn(result, Path(output_dir))
    return result
