from __future__ import annotations

import argparse
from pathlib import Path


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Keep defaults visible while allowing multi-line examples."""


def default_output_dir(input_path: str, method: str) -> Path:
    """Return the default output directory for one input-driven workflow."""

    path = Path(input_path)
    name = path.name
    stem = name[:-7] if name.endswith(".nii.gz") else path.stem
    return path.parent / f"Imt_{stem}_{method}"


def default_gene_pca_output_dir(genes_arg: str, atlas: str) -> Path:
    """Return the default gene-PCA output directory."""

    path = Path(str(genes_arg)).expanduser()
    if path.exists():
        return default_output_dir(path.as_posix(), f"gene-pca-{atlas}")
    return Path.cwd() / f"Imt_gene-pca_{atlas}"


def default_gedar_output_dir(weights_arg: str, atlas: str) -> Path:
    """Return the default GEDAR output directory."""

    path = Path(str(weights_arg)).expanduser()
    if path.exists():
        return default_output_dir(path.as_posix(), f"gedar-{atlas}")
    return Path.cwd() / f"Imt_gedar_{atlas}"


def default_gene_output_dir(gene: str, atlas: str) -> Path:
    """Return the default output directory for one single-gene query."""

    gene_name = str(gene).strip().upper() or "GENE"
    return Path.cwd() / f"Imt_gene_{gene_name}_{atlas}"


def make_shared_analysis_parent() -> argparse.ArgumentParser:
    """Build the common argument group shared by correlation and PLS runs."""

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
        choices=["default", "all", "cort+sub", "cort"],
        default="default",
        help="Which atlas regions to analyze. `default` and `cort` are cortex only; `all` and `cort+sub` include the packaged aseg add-on.",
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
        help="Gene set file or name for GSEA or ORA. Use `lake`, `pooled`, a GSEApy/Enrichr library name, or a local `.gmt` file.",
    )
    shared.add_argument(
        "--geneset-organism",
        default="Human",
        help="Organism used when --geneset is a GSEApy/Enrichr library name.",
    )
    shared.add_argument(
        "--enrichment",
        choices=["ensemble", "gsea", "ora", "none"],
        default=None,
        help=(
            "Which enrichment backend to run. `ensemble` performs phenotype-null category enrichment, "
            "`gsea` runs preranked GSEA, `ora` performs over-representation analysis, and `none` skips enrichment."
        ),
    )
    shared.add_argument(
        "--ora-p-threshold",
        type=float,
        default=None,
        help="Raw gene-level p-value threshold used when --enrichment ora is selected. Defaults to 0.05 in ORA mode.",
    )
    gsea_group = shared.add_mutually_exclusive_group()
    gsea_group.add_argument(
        "--gsea",
        dest="run_gsea",
        action="store_true",
        default=None,
        help="Legacy compatibility flag equivalent to `--enrichment gsea`.",
    )
    gsea_group.add_argument(
        "--no-gsea",
        dest="run_gsea",
        action="store_false",
        default=None,
        help="Legacy compatibility flag that disables GSEA selection. Use `--enrichment none` to skip all enrichment.",
    )
    return shared


def make_gene_pca_parent() -> argparse.ArgumentParser:
    """Build the shared parser used by the gene-PCA subcommand."""

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-g",
        "--genes",
        required=True,
        help="Gene list input. Accepts a text/TSV/CSV file or a comma-separated list of gene symbols.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output directory for tables, plots, and README.txt. Defaults to a gene-pca-specific folder.",
    )
    parser.add_argument(
        "-a",
        "--atlas",
        default="dk",
        help="Atlas preset to use. Examples: dk, schaefer-100, schaefer-200, schaefer-400, destrieux, glasser-360.",
    )
    parser.add_argument(
        "--hemisphere",
        choices=["left", "both"],
        default="left",
        help="Use the left side only, or both sides. In `both`, the right side comes from the `abagen` mirror option.",
    )
    parser.add_argument(
        "-r",
        "--regions",
        choices=["default", "all", "cort+sub", "cort"],
        default="default",
        help="Which atlas regions to analyze. `default` and `cort` are cortex only; `all` and `cort+sub` include the packaged aseg add-on.",
    )
    return parser


def make_gene_query_parent() -> argparse.ArgumentParser:
    """Build the shared parser used by the single-gene query subcommand."""

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-g",
        "--gene",
        required=True,
        help="Gene symbol to query in the selected atlas expression matrix.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output directory for tables, plots, and README.txt. Defaults to a gene-specific folder in the current directory.",
    )
    parser.add_argument(
        "-a",
        "--atlas",
        default="dk",
        help="Atlas preset to use. Examples: dk, schaefer-100, schaefer-200, schaefer-400, destrieux, glasser-360.",
    )
    parser.add_argument(
        "--hemisphere",
        choices=["left", "both"],
        default="left",
        help="Use the left side only, or both sides.",
    )
    parser.add_argument(
        "-r",
        "--regions",
        choices=["default", "all", "cort+sub", "cort"],
        default="default",
        help="Which atlas regions to analyze. `default` and `cort` are cortex only; `all` and `cort+sub` include the packaged aseg add-on.",
    )
    parser.add_argument(
        "--raw-expression",
        action="store_true",
        help="Use atlas expression values as stored instead of z-scoring the gene across regions.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=25,
        help="Maximum number of significantly positively co-expressed genes to return.",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
        help="Benjamini-Hochberg threshold used to define significant co-expression.",
    )
    return parser
