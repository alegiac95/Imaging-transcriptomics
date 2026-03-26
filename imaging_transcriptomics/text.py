from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

from .models import CorrelationResult, PLSResult


def _header_lines(result: CorrelationResult | PLSResult) -> list[str]:
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


def render_readme(result: CorrelationResult | PLSResult) -> str:
    lines = _header_lines(result)
    if isinstance(result, CorrelationResult):
        lines.extend(
            [
                "metadata.json: machine-readable run metadata, including atlas, hemisphere, and null-model settings.",
                "regional_values.tsv: region values lined up with the atlas labels.",
                "corr_genes.tsv: ranked gene table with raw p-values, BH FDR, and maxT family-wise error.",
                "plots/regional_values.png: line plot of the extracted region values.",
                "plots/corr_top_genes.png: strongest positive and negative gene hits.",
                "plots/corr_distribution.png: distribution of all gene correlations.",
            ]
        )
        if result.gsea_table is not None:
            lines.append("gsea_corr_results.tsv: GSEA results for the correlation ranking.")
            lines.append("plots/gsea_corr_dotplot.png: top GSEA terms shown as a dot plot.")
        if result.ora_tables is not None:
            lines.append("ora_corr_up.tsv / ora_corr_down.tsv: ORA results for positive and negative genes.")
            lines.append("plots/ora_corr_heatmap.png: ORA heatmap with one row for up and one row for down.")
    else:
        lines.extend(
            [
                "metadata.json: machine-readable run metadata, including atlas, hemisphere, and null-model settings.",
                "regional_values.tsv: region values lined up with the atlas labels.",
                "pls_summary.tsv: variance explained and permutation p-values for each kept PLS component.",
                "pls_component_<n>.tsv: ranked gene table for each PLS component, including weights, z-scores, BH FDR, and maxT family-wise error.",
                "plots/regional_values.png: line plot of the extracted region values.",
                "plots/pls_variance.png: variance explained by each kept PLS component.",
                "plots/pls_cumulative_variance.png: cumulative variance curve.",
                "plots/pls_component_<n>_genes.png: strongest positive and negative gene weights for each component.",
            ]
        )
        if any(component.gsea_table is not None for component in result.components):
            lines.append("gsea_pls<n>_results.tsv: GSEA results for each PLS component.")
            lines.append("plots/gsea_pls<n>_dotplot.png: top GSEA terms shown as a dot plot for each component.")
        if any(component.ora_tables is not None for component in result.components):
            lines.append("ora_pls<n>_up.tsv / ora_pls<n>_down.tsv: ORA results for positive and negative genes in each PLS component.")
            lines.append("plots/ora_pls<n>_heatmap.png: ORA heatmap with one row for up and one row for down.")
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
    return "\n".join(lines) + "\n"


def write_readme(result: CorrelationResult | PLSResult, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "README.txt"
    path.write_text(render_readme(result))
    return path
