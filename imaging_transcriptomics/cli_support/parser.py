from __future__ import annotations

import argparse

from .shared import HelpFormatter, make_gene_pca_parent, make_gene_query_parent, make_shared_analysis_parent


def _top_level_parser() -> argparse.ArgumentParser:
    """Create the root parser shared by all CLI entry points."""

    description = (
        "Run imaging transcriptomics 2.0 analyses with built-in atlas handling.\n\n"
        "Inputs can be parcel vectors, MNI152 volumetric maps, or supported surface files.\n"
        "Use `atlases` to list the atlases known to the package."
    )
    return argparse.ArgumentParser(
        description=description,
        formatter_class=HelpFormatter,
        epilog=(
            "Examples:\n"
            "  imt atlases --packaged-only\n"
            "  imt genesets --packaged-only\n"
            "  imt genesets --organism Human\n"
            "  imt corr --input map.nii.gz --space MNI152 --atlas dk --output out_dir\n"
            "  imt corr --input map.nii.gz --space MNI152 --atlas schaefer-200 \\\n"
            "    --ora-p-threshold 0.05 --geneset lake --output out_dir\n"
            "  imt pls --input map.nii.gz --space MNI152 --atlas dk --ncomp 2 --output out_dir\n"
            "  imt gene-pca --genes genes.txt --atlas dk --output out_dir\n\n"
            "  imt gene --gene RELN --atlas dk --output out_dir\n\n"
            "  imt gedar --weights twas.tsv --atlas dk --weight-column zscore --output out_dir\n\n"
            "The long command name `imagingtranscriptomics` works as well."
        ),
    )


def _add_atlases_subcommand(subparsers) -> None:
    """Register the atlas-listing subcommand."""

    atlas_parser = subparsers.add_parser(
        "atlases",
        help="List the atlases known to the package.",
        description="List the atlases known to the package, including which ones can be used right away.",
        formatter_class=HelpFormatter,
    )
    atlas_parser.add_argument(
        "--packaged-only",
        action="store_true",
        help="Show only atlases that can be used right away.",
    )


def _add_genesets_subcommand(subparsers) -> None:
    """Register the geneset-listing subcommand."""

    geneset_parser = subparsers.add_parser(
        "genesets",
        help="List packaged genesets and optional GSEApy/Enrichr libraries.",
        description=(
            "List the geneset resources accepted by --geneset. "
            "Packaged entries are always available; Enrichr libraries require gseapy and network access."
        ),
        formatter_class=HelpFormatter,
    )
    geneset_parser.add_argument(
        "--packaged-only",
        action="store_true",
        help="Show only the packaged genesets shipped with the toolbox.",
    )
    geneset_parser.add_argument(
        "--organism",
        default="Human",
        help="Organism used when listing GSEApy/Enrichr libraries.",
    )


def _add_corr_subcommand(subparsers, shared: argparse.ArgumentParser) -> None:
    """Register the spatial-correlation workflow command."""

    corr_parser = subparsers.add_parser(
        "corr",
        parents=[shared],
        help="Run the correlation workflow.",
        description=(
            "Compare one parcellated brain map with atlas gene expression, "
            "then run one enrichment backend on the ranked genes."
        ),
        formatter_class=HelpFormatter,
        epilog=(
            "Examples:\n"
            "  imt corr --input map.nii.gz --space MNI152 --atlas dk --output out_dir\n"
            "  imt corr --input map.nii.gz --space MNI152 --atlas dk \\\n"
            "    --enrichment ensemble --geneset pooled --output out_dir\n"
            "  imt corr --input map.nii.gz --space MNI152 --atlas schaefer-200 \\\n"
            "    --enrichment ora --ora-p-threshold 0.05 --geneset lake --output out_dir\n"
            "  imt corr --input map.nii.gz --space MNI152 --atlas dk \\\n"
            "    --enrichment gsea --geneset pooled --output out_dir"
        ),
    )
    corr_parser.set_defaults(method="corr")


def _add_pls_subcommand(subparsers, shared: argparse.ArgumentParser) -> None:
    """Register the PLS workflow command."""

    pls_parser = subparsers.add_parser(
        "pls",
        parents=[shared],
        help="Run the PLS workflow.",
        description=(
            "Fit a PLS model between one parcellated brain map and atlas gene expression, "
            "then run one enrichment backend on each kept component."
        ),
        formatter_class=HelpFormatter,
        epilog=(
            "Examples:\n"
            "  imt pls --input map.nii.gz --space MNI152 --atlas dk --ncomp 2 --output out_dir\n"
            "  imt pls --input map.nii.gz --space MNI152 --atlas dk --var 0.5 \\\n"
            "    --enrichment ora --ora-p-threshold 0.05 --geneset lake --output out_dir\n"
            "  imt pls --input map.nii.gz --space MNI152 --atlas dk --ncomp 2 \\\n"
            "    --enrichment ensemble --geneset GO_Biological_Process_2025 --output out_dir"
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


def _add_gene_pca_subcommand(subparsers, shared: argparse.ArgumentParser) -> None:
    """Register the gene-PCA workflow command."""

    gene_pca_parser = subparsers.add_parser(
        "gene-pca",
        parents=[shared],
        help="Run PCA on atlas expression for a selected gene list.",
        description=(
            "Filter the atlas expression matrix to the requested genes, "
            "normalize those genes across regions, and run PCA to obtain regional expression patterns."
        ),
        formatter_class=HelpFormatter,
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


def _add_gene_subcommand(subparsers, shared: argparse.ArgumentParser) -> None:
    """Register the single-gene atlas query command."""

    gene_parser = subparsers.add_parser(
        "gene",
        parents=[shared],
        help="Return one gene's atlas expression profile and top co-expressed genes.",
        description=(
            "Extract the regional expression vector for one gene from the selected atlas, "
            "plot that expression pattern, and rank the most significantly positively co-expressed genes."
        ),
        formatter_class=HelpFormatter,
        epilog=(
            "Examples:\n"
            "  imt gene --gene RELN --atlas dk --output out_dir\n"
            "  imt gene --gene GAD1 --atlas schaefer-200 --hemisphere both --top-n 50 --output out_dir\n"
            "  imt gene --gene MBP --atlas dk --raw-expression --fdr-threshold 0.01 --output out_dir"
        ),
    )
    gene_parser.set_defaults(command="gene")


def _add_gedar_subcommand(subparsers) -> None:
    """Register the GEDAR/PTRS-style weighted-expression workflow."""

    gedar_parser = subparsers.add_parser(
        "gedar",
        help="Compute a GEDAR/PTRS-style weighted average on atlas expression.",
        description=(
            "Match a weighted gene table to the selected atlas expression matrix, "
            "optionally filter genes by rank and sign, apply the packaged AHPA brain-gene filter, "
            "compute a regional GEDAR score, and optionally run GSEA or split ORA on the matched signature."
        ),
        formatter_class=HelpFormatter,
        epilog=(
            "Examples:\n"
            "  imt gedar --weights twas.tsv --atlas dk --weight-column zscore --output out_dir\n"
            "  imt gedar --weights twas.tsv --atlas schaefer-200 --rank-column fdr \\\n"
            "    --top-percent 5 --direction split --normalize-weights unit --output out_dir\n"
            "  imt gedar --weights twas.tsv --atlas dk --enrichment gsea --geneset pooled --output out_dir\n"
            "  imt gedar --weights twas.tsv --atlas dk --top-percent 10 --direction split \\\n"
            "    --enrichment ora --geneset lake --output out_dir"
        ),
    )
    gedar_parser.add_argument(
        "--weights",
        required=True,
        help="Path to a CSV/TSV table containing gene weights. Genes outside the packaged AHPA brain-gene list are removed automatically.",
    )
    gedar_parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output directory for tables, plots, and README.txt. Defaults to a GEDAR-specific folder.",
    )
    gedar_parser.add_argument(
        "-a",
        "--atlas",
        default="dk",
        help="Atlas preset to use. Examples: dk, schaefer-100, schaefer-200, schaefer-400, destrieux, glasser-360.",
    )
    gedar_parser.add_argument(
        "--hemisphere",
        choices=["left", "both"],
        default="both",
        help="Use the left side only, or both sides. In `both`, the right side comes from the `abagen` mirror option.",
    )
    gedar_parser.add_argument(
        "-r",
        "--regions",
        choices=["default", "all", "cort+sub", "cort"],
        default="default",
        help="Which atlas regions to analyze. `default` and `cort` are cortex only; `all` and `cort+sub` include the packaged aseg add-on.",
    )
    gedar_parser.add_argument(
        "--gene-column",
        default="gene",
        help="Column in the weights table containing gene symbols.",
    )
    gedar_parser.add_argument(
        "--weight-column",
        default="weight",
        help="Column in the weights table containing the gene weights used in the GEDAR average.",
    )
    gedar_parser.add_argument(
        "--rank-column",
        default=None,
        help="Optional column used to rank or threshold genes before the GEDAR average.",
    )
    gedar_parser.add_argument(
        "--rank-mode",
        choices=["ascending", "descending"],
        default="ascending",
        help="Whether smaller or larger rank values should be treated as better when selecting top genes.",
    )
    selection_group = gedar_parser.add_mutually_exclusive_group()
    selection_group.add_argument(
        "--top-percent",
        type=float,
        default=None,
        help="Keep only the top percentage of genes according to --rank-column.",
    )
    selection_group.add_argument(
        "--top-n",
        type=int,
        default=None,
        help="Keep only the top N genes according to --rank-column.",
    )
    selection_group.add_argument(
        "--p-threshold",
        type=float,
        default=None,
        help="Keep only genes passing a threshold on --rank-column.",
    )
    gedar_parser.add_argument(
        "--direction",
        choices=["combined", "up", "down", "split"],
        default="combined",
        help="Use all weighted genes together, only positive or negative weights, or return separate up and down scores in one run.",
    )
    gedar_parser.add_argument(
        "--normalize-expression",
        choices=["zscore", "none"],
        default="zscore",
        help="Whether to z-score atlas gene expression across regions before the GEDAR average.",
    )
    gedar_parser.add_argument(
        "--normalize-weights",
        choices=["none", "zscore", "unit"],
        default="none",
        help="Optional normalization applied to the selected weights before the GEDAR average.",
    )
    gedar_parser.add_argument(
        "--geneset",
        default="lake",
        help="Gene set file or name for GEDAR GSEA or ORA. Use `lake`, `pooled`, a GSEApy/Enrichr library name, or a local `.gmt` file.",
    )
    gedar_parser.add_argument(
        "--geneset-organism",
        default="Human",
        help="Organism used when --geneset is a GSEApy/Enrichr library name.",
    )
    gsea_group = gedar_parser.add_mutually_exclusive_group()
    gedar_parser.add_argument(
        "--enrichment",
        choices=["gsea", "ora", "none"],
        default=None,
        help=(
            "Optional GEDAR enrichment backend. `gsea` runs preranked GSEA on the full matched signed gene table, "
            "`ora` runs split up/down ORA on the genes selected for the GEDAR score, and `none` skips enrichment."
        ),
    )
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
        help="Legacy compatibility flag that disables GEDAR GSEA selection.",
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the full CLI parser for all supported workflows."""

    parser = _top_level_parser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    shared = make_shared_analysis_parent()
    gene_pca_shared = make_gene_pca_parent()
    gene_query_shared = make_gene_query_parent()

    _add_atlases_subcommand(subparsers)
    _add_genesets_subcommand(subparsers)
    _add_corr_subcommand(subparsers, shared)
    _add_pls_subcommand(subparsers, shared)
    _add_gene_pca_subcommand(subparsers, gene_pca_shared)
    _add_gene_subcommand(subparsers, gene_query_shared)
    _add_gedar_subcommand(subparsers)
    return parser


def parse_cmdline(argv=None):
    """Parse CLI arguments for one invocation."""

    return build_parser().parse_args(argv)
