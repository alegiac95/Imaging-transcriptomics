from __future__ import annotations

from pathlib import Path


_PACKAGED_GENESETS = {
    "lake": Path(__file__).parent / "data" / "geneset_LAKE.gmt",
    "pooled": Path(__file__).parent / "data" / "geneset_Pooled.gmt",
}


def list_packaged_genesets() -> tuple[str, ...]:
    """Return the packaged geneset names shipped with the toolbox."""

    return tuple(sorted(_PACKAGED_GENESETS))


def list_remote_genesets(organism: str = "Human") -> list[str]:
    """Return the active Enrichr library names exposed by GSEApy."""

    try:
        import gseapy
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError("gseapy is required to list remote geneset libraries.") from exc
    return sorted(gseapy.get_library_name(organism=organism))


def get_geneset(gene_set: str) -> str:
    """Return a packaged geneset path or passthrough a GSEA library name."""

    packaged = _PACKAGED_GENESETS.get(gene_set.lower())
    if packaged is not None:
        return str(packaged)
    return gene_set


def resolve_geneset_resource(gene_set: str, organism: str = "Human"):
    """Return a geneset resource suitable for both GSEA and ORA.

    The returned object is either a local GMT path or a gene-set dictionary
    downloaded through GSEApy/Enrichr.
    """

    resolved = get_geneset(gene_set)
    path = Path(resolved)
    if path.exists() and path.is_file() and path.suffix.lower() == ".gmt":
        return str(path)

    try:
        import gseapy
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "gseapy is required to resolve non-local geneset libraries. "
            f"Install gseapy or pass a local .gmt file instead of {gene_set!r}."
        ) from exc

    return gseapy.get_library(
        resolved,
        organism=organism,
        min_size=0,
        max_size=100_000,
    )
