#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from imaging_transcriptomics import atlas_table, run_corr, run_gene_pca, run_pls

from ._logging import configure_cli_logging


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Keep defaults visible while allowing multi-line examples."""


def _default_output_dir(input_path: str, method: str) -> Path:
    path = Path(input_path)
    name = path.name
    if name.endswith(".nii.gz"):
        stem = name[:-7]
    else:
        stem = path.stem
    return path.parent / f"Imt_{stem}_{method}"


def _default_gene_pca_output_dir(genes_arg: str, atlas: str) -> Path:
    path = Path(str(genes_arg)).expanduser()
    if path.exists():
        return _default_output_dir(path.as_posix(), f"gene-pca-{atlas}")
    return Path.cwd() / f"Imt_gene-pca_{atlas}"


def build_parser() -> argparse.ArgumentParser:
    description = (
        "Run imaging transcriptomics 2.0 analyses with built-in atlas handling.\n\n"
        "Inputs can be parcel vectors, MNI152 volumetric maps, or supported surface files.\n"
        "Use `atlases` to list the atlases known to the package."
    )
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=_HelpFormatter,
        epilog=(
            "Examples:\n"
            "  imt atlases --packaged-only\n"
            "  imt corr --input map.nii.gz --space MNI152 --atlas dk --output out_dir\n"
            "  imt corr --input map.nii.gz --space MNI152 --atlas schaefer-200 \\\n"
            "    --ora-p-threshold 0.05 --geneset lake --output out_dir\n"
            "  imt pls --input map.nii.gz --space MNI152 --atlas dk --ncomp 2 --output out_dir\n"
            "  imt gene-pca --genes genes.txt --atlas dk --output out_dir\n\n"
            "The long command name `imagingtranscriptomics` works as well."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    atlas_parser = subparsers.add_parser(
        "atlases",
        help="List the atlases known to the package.",
        description="List the atlases known to the package, including which ones can be used right away.",
        formatter_class=_HelpFormatter,
    )
    atlas_parser.add_argument(
        "--packaged-only",
        action="store_true",
        help="Show only atlases that can be used right away.",
    )

    gene_pca_shared = argparse.ArgumentParser(add_help=False)
    gene_pca_shared.add_argument(
        "-g",
        "--genes",
        required=True,
        help="Gene list input. Accepts a text/TSV/CSV file or a comma-separated list of gene symbols.",
    )
    gene_pca_shared.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output directory for tables, plots, and README.txt. Defaults to a gene-pca-specific folder.",
    )
    gene_pca_shared.add_argument(
        "-a",
        "--atlas",
        default="dk",
        help="Atlas preset to use. Examples: dk, schaefer-100, schaefer-200, schaefer-400, destrieux, glasser-360.",
    )
    gene_pca_shared.add_argument(
        "--hemisphere",
        choices=["left", "both"],
        default="left",
        help="Use the left side only, or both sides. In `both`, the right side comes from the `abagen` mirror option.",
    )
    gene_pca_shared.add_argument(
        "-r",
        "--regions",
        choices=["all", "cort+sub", "cort"],
        default="all",
        help="Which atlas regions to analyze: all regions, cortex plus subcortex, or cortex only.",
    )

    shared = argparse.ArgumentParser(add_help=False)
    shared.add_argument(
        "-i",
        "--input",
        required=True,
        help=(
            "Input data. Accepts a parcel vector text file, a volumetric NIfTI map, or a supported surface file. "
            "Raw subject-space anatomical scans should be registered to MNI152 first."
        ),
    )
    shared.add_argument(
        "--input-rh",
        default=None,
        help="Right-hemisphere surface file to pair with a left-hemisphere surface input.",
    )
    shared.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output directory for tables, plots, and README.txt. Defaults to a method-specific folder next to the input.",
    )
    shared.add_argument(
        "-a",
        "--atlas",
        default="dk",
        help="Atlas preset to use. Examples: dk, schaefer-100, schaefer-200, schaefer-400, destrieux, glasser-360.",
    )
    shared.add_argument(
        "--hemisphere",
        choices=["left", "both"],
        default="left",
        help="Use the left side only, or both sides. In `both`, the right side comes from the `abagen` mirror option.",
    )
    shared.add_argument(
        "-r",
        "--regions",
        choices=["all", "cort+sub", "cort"],
        default="all",
        help="Which atlas regions to analyze: all regions, cortex plus subcortex, or cortex only.",
    )
    shared.add_argument(
        "--space",
        default=None,
        help="Source space for the input map, for example MNI152, fsaverage, fsLR, or CIVET.",
    )
    shared.add_argument(
        "-p",
        "--permutations",
        type=int,
        default=1000,
        help="Number of permutations or null samples to use.",
    )
    shared.add_argument(
        "--seed",
        type=int,
        default=1234,
        help="Random seed for permutation and null-model generation.",
    )
    shared.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of worker threads to use for parallel work such as PLS permutation fits.",
    )
    shared.add_argument(
        "--null-method",
        choices=["auto", "vasa", "alexander_bloch", "moran", "random"],
        default="auto",
        help=(
            "How cortical nulls should be generated. "
            "'auto' tries Vasa spins first and falls back to random shuffling within each hemisphere."
        ),
    )
    shared.add_argument(
        "--geneset",
        default="lake",
        help="Gene set file or name for GSEA or ORA. Use `lake`, `pooled`, or a local `.gmt` file.",
    )
    shared.add_argument(
        "--ora-p-threshold",
        type=float,
        default=None,
        help="Run ORA on genes with raw p-values at or below this threshold, separately for positive and negative scores.",
    )
    gsea_group = shared.add_mutually_exclusive_group()
    gsea_group.add_argument(
        "--gsea",
        dest="run_gsea",
        action="store_true",
        default=None,
        help="Force GSEA on. Useful when you also request ORA and want both analyses.",
    )
    gsea_group.add_argument(
        "--no-gsea",
        dest="run_gsea",
        action="store_false",
        default=None,
        help="Skip GSEA and only write the regional, gene, and ORA outputs.",
    )

    corr_parser = subparsers.add_parser(
        "corr",
        parents=[shared],
        help="Run the correlation workflow.",
        description=(
            "Compare one parcellated brain map with atlas gene expression, "
            "then optionally run GSEA and/or ORA on the ranked genes."
        ),
        formatter_class=_HelpFormatter,
        epilog=(
            "Examples:\n"
            "  imt corr --input map.nii.gz --space MNI152 --atlas dk --output out_dir\n"
            "  imt corr --input map.nii.gz --space MNI152 --atlas schaefer-200 \\\n"
            "    --ora-p-threshold 0.05 --geneset lake --output out_dir\n"
            "  imt corr --input map.nii.gz --space MNI152 --atlas dk \\\n"
            "    --ora-p-threshold 0.05 --gsea --geneset pooled --output out_dir"
        ),
    )
    corr_parser.set_defaults(method="corr")

    pls_parser = subparsers.add_parser(
        "pls",
        parents=[shared],
        help="Run the PLS workflow.",
        description=(
            "Fit a PLS model between one parcellated brain map and atlas gene expression, "
            "then optionally run GSEA and/or ORA on each kept component."
        ),
        formatter_class=_HelpFormatter,
        epilog=(
            "Examples:\n"
            "  imt pls --input map.nii.gz --space MNI152 --atlas dk --ncomp 2 --output out_dir\n"
            "  imt pls --input map.nii.gz --space MNI152 --atlas dk --var 0.5 \\\n"
            "    --ora-p-threshold 0.05 --geneset lake --output out_dir"
        ),
    )
    pls_group = pls_parser.add_mutually_exclusive_group(required=True)
    pls_group.add_argument(
        "--ncomp",
        type=int,
        default=None,
        help="Number of PLS components to retain.",
    )
    pls_group.add_argument(
        "--var",
        type=float,
        default=None,
        help="Cumulative variance target used to choose how many PLS components to keep.",
    )
    pls_parser.set_defaults(method="pls")

    gene_pca_parser = subparsers.add_parser(
        "gene-pca",
        parents=[gene_pca_shared],
        help="Run PCA on atlas expression for a selected gene list.",
        description=(
            "Filter the atlas expression matrix to the requested genes, "
            "normalize those genes across regions, and run PCA to obtain regional expression patterns."
        ),
        formatter_class=_HelpFormatter,
        epilog=(
            "Examples:\n"
            "  imt gene-pca --genes genes.txt --atlas dk --output out_dir\n"
            "  imt gene-pca --genes RELN,GAD1,SLC1A2 --atlas schaefer-100 --hemisphere both --ncomp 2 --output out_dir"
        ),
    )
    gene_pca_parser.add_argument(
        "--ncomp",
        type=int,
        default=3,
        help="Maximum number of PCA components to retain.",
    )
    return parser


def parse_cmdline(argv=None):
    return build_parser().parse_args(argv)


def _resolve_run_gsea(parsed) -> bool:
    if parsed.run_gsea is not None:
        return bool(parsed.run_gsea)
    return parsed.ora_p_threshold is None


def main():
    parsed = parse_cmdline()
    configure_cli_logging()
    if parsed.command == "atlases":
        print(atlas_table(packaged_only=parsed.packaged_only).to_string(index=False))
        return

    if parsed.command == "gene-pca":
        output_dir = Path(parsed.output) if parsed.output else _default_gene_pca_output_dir(parsed.genes, parsed.atlas)
        run_gene_pca(
            parsed.genes,
            atlas=parsed.atlas,
            hemisphere=parsed.hemisphere,
            regions=parsed.regions,
            n_components=parsed.ncomp,
            output_dir=output_dir,
        )
        return

    output_dir = Path(parsed.output) if parsed.output else _default_output_dir(parsed.input, parsed.method)
    run_gsea = _resolve_run_gsea(parsed)

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
        ora_p_threshold=parsed.ora_p_threshold,
        seed=parsed.seed,
        n_jobs=parsed.jobs,
    )


if __name__ == "__main__":
    main()
