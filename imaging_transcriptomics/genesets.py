from __future__ import annotations

from pathlib import Path


def get_geneset(gene_set: str) -> str:
    """Return a packaged geneset path or passthrough a GSEA library name."""

    if gene_set.lower() == "lake":
        return str(Path(__file__).parent / "data" / "geneset_LAKE.gmt")
    if gene_set.lower() == "pooled":
        return str(Path(__file__).parent / "data" / "geneset_Pooled.gmt")
    return gene_set
