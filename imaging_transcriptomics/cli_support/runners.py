from __future__ import annotations

from pathlib import Path

from imaging_transcriptomics import atlas_table, run_corr, run_gedar, run_gene_pca, run_pls

from ..genesets import list_packaged_genesets, list_remote_genesets
from .shared import default_gene_pca_output_dir, default_gedar_output_dir, default_output_dir


def resolve_run_gsea(parsed) -> bool:
    """Resolve the default GSEA behavior from CLI flags."""

    if parsed.run_gsea is not None:
        return bool(parsed.run_gsea)
    return parsed.ora_p_threshold is None


def run_atlases_command(parsed) -> None:
    """Handle the `imt atlases` subcommand."""

    print(atlas_table(packaged_only=parsed.packaged_only).to_string(index=False))


def run_genesets_command(parsed) -> None:
    """Handle the `imt genesets` subcommand."""

    packaged = list_packaged_genesets()
    if parsed.packaged_only:
        print("Packaged genesets")
        for name in packaged:
            print(name)
        return

    print("Packaged genesets")
    for name in packaged:
        print(name)
    print("")
    print(f"GSEApy/Enrichr genesets ({parsed.organism})")
    try:
        remote = list_remote_genesets(parsed.organism)
    except Exception as exc:
        print(f"Unable to list remote genesets: {exc}")
        return
    for name in remote:
        print(name)


def run_gene_pca_command(parsed) -> None:
    """Handle the `imt gene-pca` subcommand."""

    output_dir = Path(parsed.output) if parsed.output else default_gene_pca_output_dir(parsed.genes, parsed.atlas)
    run_gene_pca(
        parsed.genes,
        atlas=parsed.atlas,
        hemisphere=parsed.hemisphere,
        regions=parsed.regions,
        n_components=parsed.ncomp,
        output_dir=output_dir,
    )


def run_gedar_command(parsed) -> None:
    """Handle the `imt gedar` subcommand."""

    output_dir = Path(parsed.output) if parsed.output else default_gedar_output_dir(parsed.weights, parsed.atlas)
    run_gedar(
        parsed.weights,
        atlas=parsed.atlas,
        hemisphere=parsed.hemisphere,
        regions=parsed.regions,
        gene_column=parsed.gene_column,
        weight_column=parsed.weight_column,
        rank_column=parsed.rank_column,
        rank_mode=parsed.rank_mode,
        top_percent=parsed.top_percent,
        top_n=parsed.top_n,
        p_threshold=parsed.p_threshold,
        direction=parsed.direction,
        normalize_expression=parsed.normalize_expression,
        normalize_weights=parsed.normalize_weights,
        output_dir=output_dir,
    )


def run_analysis_command(parsed) -> None:
    """Handle the `imt corr` and `imt pls` workflow subcommands."""

    output_dir = Path(parsed.output) if parsed.output else default_output_dir(parsed.input, parsed.method)
    run_gsea = resolve_run_gsea(parsed)

    if parsed.method == "corr":
        run_corr(
            parsed.input,
            atlas=parsed.atlas,
            hemisphere=parsed.hemisphere,
            regions=parsed.regions,
            source_space=parsed.space,
            input_rh=parsed.input_rh,
            n_permutations=parsed.permutations,
            null_method=parsed.null_method,
            output_dir=output_dir,
            run_gsea=run_gsea,
            gene_set=parsed.geneset,
            geneset_organism=parsed.geneset_organism,
            ora_p_threshold=parsed.ora_p_threshold,
            seed=parsed.seed,
            n_jobs=parsed.jobs,
        )
        return

    run_pls(
        parsed.input,
        atlas=parsed.atlas,
        hemisphere=parsed.hemisphere,
        regions=parsed.regions,
        source_space=parsed.space,
        input_rh=parsed.input_rh,
        n_components=parsed.ncomp,
        var=parsed.var,
        n_permutations=parsed.permutations,
        null_method=parsed.null_method,
        output_dir=output_dir,
        run_gsea=run_gsea,
        gene_set=parsed.geneset,
        geneset_organism=parsed.geneset_organism,
        ora_p_threshold=parsed.ora_p_threshold,
        seed=parsed.seed,
        n_jobs=parsed.jobs,
    )


def dispatch_command(parsed) -> None:
    """Dispatch one parsed CLI namespace to the matching workflow."""

    if parsed.command == "atlases":
        run_atlases_command(parsed)
        return
    if parsed.command == "genesets":
        run_genesets_command(parsed)
        return
    if parsed.command == "gene-pca":
        run_gene_pca_command(parsed)
        return
    if parsed.command == "gedar":
        run_gedar_command(parsed)
        return
    run_analysis_command(parsed)
