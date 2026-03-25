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
        f"Spatial null model: {meta.null_method}",
    ]
    if meta.n_components is not None:
        lines.append(f"PLS components: {meta.n_components}")
    if meta.geneset is not None:
        lines.append(f"GSEA geneset: {meta.geneset}")
    lines.extend(["", "Files", "-----"])
    return lines


def render_readme(result: CorrelationResult | PLSResult) -> str:
    lines = _header_lines(result)
    if isinstance(result, CorrelationResult):
        lines.extend(
            [
                "metadata.json: machine-readable run metadata, including atlas, hemisphere, and null-model settings.",
                "regional_values.tsv: parcellated scan values aligned to the atlas labels.",
                "corr_genes.tsv: ranked correlation table with raw and FDR-corrected p-values.",
                "plots/regional_values.png: profile of the extracted regional scan values.",
                "plots/corr_top_genes.png: strongest positive and negative gene associations.",
                "plots/corr_distribution.png: distribution of all gene-wise correlations.",
            ]
        )
        if result.gsea_table is not None:
            lines.append("gsea_corr_results.tsv: optional gene-set enrichment results.")
            lines.append("plots/gsea_corr_dotplot.png: top enriched pathways shown as a GSEA dot plot.")
    else:
        lines.extend(
            [
                "metadata.json: machine-readable run metadata, including atlas, hemisphere, and null-model settings.",
                "regional_values.tsv: parcellated scan values aligned to the atlas labels.",
                "pls_summary.tsv: per-component explained variance and permutation p-values.",
                "pls_component_<n>.tsv: ranked gene table for each PLS component, including weights and z-scores.",
                "plots/regional_values.png: profile of the extracted regional scan values.",
                "plots/pls_variance.png: variance explained by each retained PLS component.",
                "plots/pls_cumulative_variance.png: cumulative explained variance curve.",
                "plots/pls_component_<n>_genes.png: strongest positive and negative gene weights per component.",
            ]
        )
        if any(component.gsea_table is not None for component in result.components):
            lines.append("gsea_pls<n>_results.tsv: optional gene-set enrichment results per component.")
            lines.append("plots/gsea_pls<n>_dotplot.png: top enriched pathways shown as a GSEA dot plot for each component.")
    lines.extend(
        [
            "",
            "Notes",
            "-----",
            "This v2 branch keeps the abagen-derived expression matrices and adds hemisphere-aware atlas selection.",
            "Packaged atlases now ship ready-to-run assets for dk, schaefer-100, schaefer-200, schaefer-400, destrieux, and glasser-360.",
            "When neuromaps is installed, non-MNI inputs and surface inputs can be resampled/parcellated through the new scan extraction helpers.",
        ]
    )
    return "\n".join(lines) + "\n"


def write_readme(result: CorrelationResult | PLSResult, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "README.txt"
    path.write_text(render_readme(result))
    return path
