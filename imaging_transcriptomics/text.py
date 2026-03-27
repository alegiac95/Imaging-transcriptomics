from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

from .models import CorrelationResult, GenePCAResult, PLSResult


def _analysis_header_lines(result: CorrelationResult | PLSResult) -> list[str]:
    meta = result.metadata
    lines = [
        "Imaging Transcriptomics 2.0 Run",
        "==============================",
        "",
        f"Generated: {datetime.now(timezone.utc).astimezone().isoformat(timespec='seconds')}",
        f"Method: {meta.method}",
        f"Atlas: {meta.atlas_label} ({meta.atlas_id})",
        f"Hemisphere mode: {meta.hemisphere}",
        f"Region scope: {meta.regions}",
        f"Input source: {meta.source}",
        f"Input kind: {meta.source_kind}",
        f"Input space: {meta.source_space or 'unspecified'}",
        f"Permutations: {meta.n_permutations}",
        f"Null method: {meta.null_method}",
    ]
    if meta.n_components is not None:
        lines.append(f"PLS components kept: {meta.n_components}")
    if meta.geneset is not None:
        lines.append(f"Gene set file or name: {meta.geneset}")
    if meta.ora_p_threshold is not None:
        lines.append(f"ORA p-value threshold: {meta.ora_p_threshold}")
    lines.extend(["", "Files", "-----"])
    return lines


def _gene_pca_header_lines(result: GenePCAResult) -> list[str]:
    return [
        "Imaging Transcriptomics 2.0 Gene PCA",
        "===================================",
        "",
        f"Generated: {datetime.now(timezone.utc).astimezone().isoformat(timespec='seconds')}",
        "Method: gene-pca",
        f"Atlas: {result.atlas_label} ({result.atlas_id})",
        f"Hemisphere mode: {result.hemisphere}",
        f"Region scope: {result.regions}",
        f"Requested genes: {len(result.requested_genes)}",
        f"Genes used in PCA: {len(result.matched_genes)}",
        f"Missing genes: {len(result.missing_genes)}",
        "",
        "Files",
        "-----",
    ]


def render_readme(result: CorrelationResult | PLSResult | GenePCAResult, *, plots_available: bool = True) -> str:
    if isinstance(result, GenePCAResult):
        lines = _gene_pca_header_lines(result)
        lines.extend(
            [
                "metadata.json: machine-readable metadata for the gene PCA run.",
                "gene_pca_scores.tsv: regional PCA scores for each retained component.",
                "gene_pca_loadings.tsv: per-gene PCA loadings for each retained component.",
                "gene_pca_variance.tsv: variance explained and cumulative variance by component.",
                "matched_genes.txt: genes from the input list that were found in the atlas expression matrix.",
                "missing_genes.txt: input genes that were not found in the atlas expression matrix.",
            ]
        )
        if plots_available:
            lines.extend(
                [
                    "plots/gene_pca_variance.png: variance explained by each PCA component.",
                    "plots/gene_pca_pc<n>_regions.png: regional score profile for each plotted component.",
                    "plots/gene_pca_pc<n>_loadings.png: strongest positive and negative gene loadings for each plotted component.",
                ]
            )
        lines.extend(
            [
                "",
                "Notes",
                "-----",
                "Gene expression values are standardized across regions before PCA.",
                "PCA is run on the selected atlas expression matrix after filtering to the requested genes.",
            ]
        )
        if not plots_available:
            lines.append("Plot PNGs were skipped because the optional plotting dependencies are not installed.")
        return "\n".join(lines) + "\n"

    lines = _analysis_header_lines(result)
    if isinstance(result, CorrelationResult):
        lines.extend(
            [
                "metadata.json: machine-readable run metadata, including atlas, hemisphere, and null-model settings.",
                "regional_values.tsv: region values lined up with the atlas labels.",
                "corr_genes.tsv: ranked gene table with raw p-values, BH FDR, and maxT family-wise error.",
            ]
        )
        if plots_available:
            lines.extend(
                [
                    "plots/regional_values.png: line plot of the extracted region values.",
                    "plots/corr_top_genes.png: strongest positive and negative gene hits.",
                    "plots/corr_distribution.png: distribution of all gene correlations.",
                ]
            )
        if result.gsea_table is not None:
            lines.append("gsea_corr_results.tsv: GSEA results for the correlation ranking.")
            if plots_available:
                lines.append("plots/gsea_corr_dotplot.png: top GSEA terms shown as a dot plot.")
        if result.ora_tables is not None:
            lines.append("ora_corr_up.tsv / ora_corr_down.tsv: ORA results for positive and negative genes, including odds ratios and 95% confidence intervals.")
            if plots_available:
                lines.append("plots/ora_corr_heatmap.png: ORA heatmap with one row for up and one row for down, annotated with odds ratios and significance stars.")
    else:
        lines.extend(
            [
                "metadata.json: machine-readable run metadata, including atlas, hemisphere, and null-model settings.",
                "regional_values.tsv: region values lined up with the atlas labels.",
                "pls_summary.tsv: variance explained and permutation p-values for each kept PLS component.",
                "pls_component_<n>.tsv: ranked gene table for each PLS component, including weights, z-scores, BH FDR, and maxT family-wise error.",
            ]
        )
        if plots_available:
            lines.extend(
                [
                    "plots/regional_values.png: line plot of the extracted region values.",
                    "plots/pls_variance.png: variance explained by each kept PLS component.",
                    "plots/pls_cumulative_variance.png: cumulative variance curve.",
                    "plots/pls_component_<n>_genes.png: strongest positive and negative gene weights for each component.",
                ]
            )
        if any(component.gsea_table is not None for component in result.components):
            lines.append("gsea_pls<n>_results.tsv: GSEA results for each PLS component.")
            if plots_available:
                lines.append("plots/gsea_pls<n>_dotplot.png: top GSEA terms shown as a dot plot for each component.")
        if any(component.ora_tables is not None for component in result.components):
            lines.append("ora_pls<n>_up.tsv / ora_pls<n>_down.tsv: ORA results for positive and negative genes in each PLS component, including odds ratios and 95% confidence intervals.")
            if plots_available:
                lines.append("plots/ora_pls<n>_heatmap.png: ORA heatmap with one row for up and one row for down, annotated with odds ratios and significance stars.")
    lines.extend(
        [
            "",
            "Notes",
            "-----",
            "This version keeps the `abagen`-derived expression matrices and supports left-only or both-hemisphere atlas data.",
            "Included atlases are dk, schaefer-100, schaefer-200, schaefer-400, destrieux, and glasser-360.",
            "When `neuromaps` is installed, non-MNI and surface inputs can be resampled before analysis.",
        ]
    )
    if not plots_available:
        lines.append("Plot PNGs were skipped because the optional plotting dependencies are not installed.")
    return "\n".join(lines) + "\n"


def write_readme(
    result: CorrelationResult | PLSResult | GenePCAResult,
    output_dir: Path,
    *,
    plots_available: bool = True,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "README.txt"
    path.write_text(render_readme(result, plots_available=plots_available))
    return path
