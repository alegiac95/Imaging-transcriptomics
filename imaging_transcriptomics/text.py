from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

from .models import CorrelationResult, GEDARResult, GenePCAResult, PLSResult


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
    if meta.geneset_organism is not None:
        lines.append(f"Gene set organism: {meta.geneset_organism}")
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
        f"Filtered by brain gene list: {len(result.brain_filtered_genes)}",
        f"Missing genes: {len(result.missing_genes)}",
        "",
        "Files",
        "-----",
    ]


def _gedar_header_lines(result: GEDARResult) -> list[str]:
    if result.direction == "split":
        selected_line = (
            f"Selected genes: up={int(result.gene_table['selected_up'].astype(bool).sum())}, "
            f"down={int(result.gene_table['selected_down'].astype(bool).sum())}"
        )
    else:
        selected_line = f"Selected genes: {int(result.gene_table['selected'].astype(bool).sum())}"
    return [
        "Imaging Transcriptomics 2.0 GEDAR",
        "================================",
        "",
        f"Generated: {datetime.now(timezone.utc).astimezone().isoformat(timespec='seconds')}",
        "Method: gedar",
        f"Atlas: {result.atlas_label} ({result.atlas_id})",
        f"Hemisphere mode: {result.hemisphere}",
        f"Region scope: {result.regions}",
        f"Weight source: {result.weights_source}",
        f"Direction filter: {result.direction}",
        f"Expression normalization: {result.normalize_expression}",
        f"Weight normalization: {result.normalize_weights}",
        "Brain gene filter: packaged AHPA_mrna_brain.tsv",
        f"Requested genes: {len(result.requested_genes)}",
        f"Matched genes: {len(result.matched_genes)}",
        f"Missing genes: {len(result.missing_genes)}",
        selected_line,
        f"Excluded input rows: {int(result.excluded_table.shape[0])}",
        "",
        "Files",
        "-----",
    ]


def render_readme(result: CorrelationResult | PLSResult | GenePCAResult | GEDARResult, *, plots_available: bool = True) -> str:
    if isinstance(result, GenePCAResult):
        lines = _gene_pca_header_lines(result)
        lines.extend(
            [
                "metadata.json: machine-readable metadata for the gene PCA run.",
                "gene_pca_scores.tsv: regional PCA scores for each retained component.",
                "gene_pca_loadings.tsv: per-gene PCA loadings for each retained component.",
                "gene_pca_variance.tsv: variance explained and cumulative variance by component.",
                "matched_genes.txt: genes from the input list that were found in the atlas expression matrix.",
                "brain_filtered_genes.txt: input genes removed by the packaged AHPA brain-gene filter before atlas matching.",
                "missing_genes.txt: input genes that were not found in the atlas expression matrix.",
            ]
        )
        if plots_available:
            lines.extend(
                [
                    "plots/gene_pca_variance.png: variance explained by each PCA component.",
                    "plots/gene_pca_pc<n>_brain.png: atlas brain map of the regional PCA score for each plotted component.",
                    "plots/gene_pca_pc<n>_cortex.png: publication-oriented cortical surface map of the regional PCA score when surface atlas geometry is available.",
                    "plots/gene_pca_pc<n>_cortex_brainspace.png: optional BrainSpace-rendered cortical comparison plot when the BrainSpace backend is available.",
                    "plots/gene_pca_pc<n>_regions.png: regional score profile for each plotted component.",
                    "plots/gene_pca_pc<n>_loadings.png: strongest positive and negative gene loadings for each plotted component.",
                ]
            )
        lines.extend(
            [
                "",
                "Notes",
                "-----",
                "The packaged AHPA brain-gene list is applied before atlas matching, and any removed genes are written to `brain_filtered_genes.txt`.",
                "Gene expression values are standardized across regions before PCA.",
                "PCA is run on the selected atlas expression matrix after filtering to the requested genes.",
            ]
        )
        if not plots_available:
            lines.append("Plot PNGs were skipped because the optional plotting dependencies are not installed.")
        return "\n".join(lines) + "\n"

    if isinstance(result, GEDARResult):
        lines = _gedar_header_lines(result)
        lines.extend(
            [
                "metadata.json: machine-readable metadata for the GEDAR run.",
                "gedar_scores.tsv: regional weighted-average expression scores. Split runs contain separate up/down score columns.",
                "gedar_genes.tsv: matched gene weights, rank values, and selection flags used in the GEDAR average.",
                "gedar_excluded.tsv: input rows excluded before atlas matching, including duplicate symbols and invalid values.",
                "matched_genes.txt: input genes that were found in the atlas expression matrix.",
                "missing_genes.txt: input genes that were not found in the atlas expression matrix.",
            ]
        )
        if plots_available and result.direction != "split":
            lines.extend(
                [
                    "plots/gedar_scores.png: regional GEDAR score profile across the selected atlas regions.",
                    "plots/gedar_brain.png: atlas brain map of the z-scored GEDAR regional score.",
                    "plots/gedar_cortex.png: publication-oriented cortical surface map of the z-scored GEDAR regional score when surface atlas geometry is available.",
                    "plots/gedar_cortex_brainspace.png: optional BrainSpace-rendered cortical comparison plot when the BrainSpace backend is available.",
                    "plots/gedar_weights.png: strongest positive and negative gene weights used in the GEDAR average.",
                ]
            )
        if plots_available and result.direction == "split":
            lines.extend(
                [
                    "plots/gedar_up_scores.png / plots/gedar_down_scores.png: regional GEDAR score profiles for the separate up and down runs.",
                    "plots/gedar_up_brain.png / plots/gedar_down_brain.png: atlas brain maps of the z-scored up and down GEDAR scores.",
                    "plots/gedar_up_cortex.png / plots/gedar_down_cortex.png: publication-oriented cortical surface maps of the z-scored up and down GEDAR scores when surface atlas geometry is available.",
                    "plots/gedar_up_cortex_brainspace.png / plots/gedar_down_cortex_brainspace.png: optional BrainSpace-rendered cortical comparison plots when the BrainSpace backend is available.",
                    "plots/gedar_up_weights.png / plots/gedar_down_weights.png: strongest gene weights used in the separate up and down GEDAR averages.",
                ]
            )
        lines.extend(
            [
                "",
                "Notes",
                "-----",
                "GEDAR is implemented here as an atlas-agnostic weighted average of selected gene weights on the packaged regional expression matrix.",
                "Rows with blank genes, non-finite weights, non-finite rank values, duplicated gene symbols, or genes outside the packaged AHPA brain-gene filter are removed automatically and recorded in `gedar_excluded.tsv`.",
                "In `combined` mode the signed weights are averaged as in the original PTRS script. In `up` and `down` modes the selected weights are converted to absolute values before averaging, again matching the original PTRS workflow.",
                "The output `score` is the regional weighted average across genes, while `score_z` is standardized across the selected regions. Split runs write `score_up`, `score_up_z`, `score_down`, and `score_down_z` instead.",
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
                "corr_genes.tsv: ranked gene table with score, p, FDR, and maxT.",
            ]
        )
        if plots_available:
            lines.extend(
                [
                    "plots/regional_values_brain.png: atlas brain map of the extracted regional values.",
                    "plots/regional_values_cortex.png: publication-oriented cortical surface map of z-scored regional values when surface atlas geometry is available.",
                    "plots/regional_values_cortex_brainspace.png: optional BrainSpace-rendered cortical comparison plot when the BrainSpace backend is available.",
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
                lines.append("plots/ora_corr_heatmap.png: ORA heatmap with one row for up and one row for down, colored by enrichment significance and annotated with odds ratios and significance stars.")
    else:
        lines.extend(
            [
                "metadata.json: machine-readable run metadata, including atlas, hemisphere, and null-model settings.",
                "regional_values.tsv: region values lined up with the atlas labels.",
                "pls_summary.tsv: variance explained, cumulative variance, and component p values.",
                "pls_component_<n>.tsv: ranked gene table for each PLS component, including weight, zscore, p, FDR, and maxT.",
            ]
        )
        if plots_available:
            lines.extend(
                [
                    "plots/regional_values_brain.png: atlas brain map of the extracted regional values.",
                    "plots/regional_values_cortex.png: publication-oriented cortical surface map of z-scored regional values when surface atlas geometry is available.",
                    "plots/regional_values_cortex_brainspace.png: optional BrainSpace-rendered cortical comparison plot when the BrainSpace backend is available.",
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
                lines.append("plots/ora_pls<n>_heatmap.png: ORA heatmap with one row for up and one row for down, colored by enrichment significance and annotated with odds ratios and significance stars.")
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
        lines.append(
            "Plot PNGs were skipped because Matplotlib is unavailable in the "
            "current environment."
        )
    return "\n".join(lines) + "\n"


def write_readme(
    result: CorrelationResult | PLSResult | GenePCAResult | GEDARResult,
    output_dir: Path,
    *,
    plots_available: bool = True,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "README.txt"
    path.write_text(render_readme(result, plots_available=plots_available))
    return path
