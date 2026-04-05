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


def parse_gmt(path: str | Path) -> dict[str, tuple[str, ...]]:
    """Parse a GMT file into a mapping of term name to gene symbols."""

    gmt_path = Path(path)
    genesets: dict[str, tuple[str, ...]] = {}
    with gmt_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            term = fields[0]
            genes = tuple(gene for gene in fields[2:] if gene)
            genesets[term] = genes
    return genesets


def as_geneset_mapping(resource) -> dict[str, tuple[str, ...]]:
    """Normalize a geneset resource into a plain term-to-genes mapping."""

    if isinstance(resource, (str, Path)):
        return parse_gmt(resource)
    if isinstance(resource, dict):
        mapping: dict[str, tuple[str, ...]] = {}
        for term, genes in resource.items():
            if isinstance(genes, str):
                mapping[str(term)] = (genes,)
            else:
                mapping[str(term)] = tuple(str(gene) for gene in genes)
        return mapping
    raise TypeError(
        "Geneset resource must be a GMT path or a mapping of term names to iterables of gene symbols."
    )
