"""Public facade for the GEDAR weighted-expression workflow."""

from __future__ import annotations

from pathlib import Path

from .gene_expression import load_brain_gene_symbols, select_atlas_data
from .serialization import write_result_bundle
from .workflows.gedar import run_gedar as _run_gedar


def run_gedar(
    weights,
    *,
    atlas: str = "dk",
    hemisphere: str = "both",
    regions: str = "all",
    gene_column: str = "gene",
    weight_column: str = "weight",
    rank_column: str | None = None,
    rank_mode: str = "ascending",
    top_percent: float | None = None,
    top_n: int | None = None,
    p_threshold: float | None = None,
    direction: str = "combined",
    normalize_expression: str = "zscore",
    normalize_weights: str = "none",
    output_dir: str | Path | None = None,
):
    """Compute a PTRS-style weighted gene-expression score on a packaged atlas."""

    return _run_gedar(
        weights,
        atlas=atlas,
        hemisphere=hemisphere,
        regions=regions,
        gene_column=gene_column,
        weight_column=weight_column,
        rank_column=rank_column,
        rank_mode=rank_mode,
        top_percent=top_percent,
        top_n=top_n,
        p_threshold=p_threshold,
        direction=direction,
        normalize_expression=normalize_expression,
        normalize_weights=normalize_weights,
        output_dir=output_dir,
        select_atlas_data_fn=select_atlas_data,
        load_brain_gene_symbols_fn=load_brain_gene_symbols,
        write_result_bundle_fn=write_result_bundle,
    )
