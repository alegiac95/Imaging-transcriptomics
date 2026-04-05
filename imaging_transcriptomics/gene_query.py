"""Public facade for the single-gene atlas query workflow."""

from __future__ import annotations

from pathlib import Path

from .gene_expression import select_atlas_data
from .serialization import write_result_bundle
from .workflows.genes import run_gene_query as _run_gene_query


def run_gene(
    gene: str,
    *,
    atlas: str = "dk",
    hemisphere: str = "left",
    regions: str = "default",
    zscore_expression: bool = True,
    top_n: int = 25,
    fdr_threshold: float = 0.05,
    output_dir: str | Path | None = None,
):
    """Return one gene's regional expression vector and top co-expressed genes."""

    return _run_gene_query(
        gene,
        atlas=atlas,
        hemisphere=hemisphere,
        regions=regions,
        zscore_expression=zscore_expression,
        top_n=top_n,
        fdr_threshold=fdr_threshold,
        output_dir=output_dir,
        select_atlas_data_fn=select_atlas_data,
        write_result_bundle_fn=write_result_bundle,
    )
