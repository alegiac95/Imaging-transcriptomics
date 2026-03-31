"""Public facade for the gene-PCA workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from .gene_expression import select_atlas_data
from .serialization import write_result_bundle
from .workflows.gene_pca import run_gene_pca as _run_gene_pca


def run_gene_pca(
    genes: str | Path | Iterable[str],
    *,
    atlas: str = "dk",
    hemisphere: str = "left",
    regions: str = "default",
    n_components: int = 3,
    output_dir: str | Path | None = None,
):
    """Run PCA on atlas expression after filtering to a selected gene list."""

    return _run_gene_pca(
        genes,
        atlas=atlas,
        hemisphere=hemisphere,
        regions=regions,
        n_components=n_components,
        output_dir=output_dir,
        select_atlas_data_fn=select_atlas_data,
        write_result_bundle_fn=write_result_bundle,
    )
