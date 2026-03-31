from __future__ import annotations

from pathlib import Path

import numpy as np

from ..models import CorrelationResult, GEDARResult, GenePCAResult, PLSResult
from .brain import (
    plot_brain_volume_map,
    plot_cortical_surface_map,
    plot_cortical_surface_map_brainspace,
)
from .common import zscore_for_plot
from .enrichment import plot_gsea_dotplot, plot_ora_heatmap
from .gene_tables import (
    plot_correlation_distribution,
    plot_correlation_ranking,
    plot_gedar_regional_scores,
    plot_gedar_weights,
    plot_gene_pca_loadings,
    plot_gene_pca_regional_component,
    plot_gene_pca_variance,
    plot_pls_component,
    plot_pls_variance,
)


def _append_cortical_plots(
    paths: list[Path],
    *,
    table,
    atlas_id: str,
    value_column: str,
    title: str,
    output_path: Path,
) -> None:
    """Write the default cortical plot plus the optional BrainSpace variant."""

    cortex_path = plot_cortical_surface_map(
        table,
        atlas_id=atlas_id,
        value_column=value_column,
        title=title,
        output_path=output_path,
    )
    if cortex_path is not None:
        paths.append(cortex_path)

    brainspace_path = plot_cortical_surface_map_brainspace(
        table,
        atlas_id=atlas_id,
        value_column=value_column,
        title=title,
        output_path=output_path.with_name(f"{output_path.stem}_brainspace{output_path.suffix}"),
    )
    if brainspace_path is not None:
        paths.append(brainspace_path)


def save_result_plots(result: CorrelationResult | PLSResult | GenePCAResult | GEDARResult, output_dir: Path) -> list[Path]:
    """Write the standard plot bundle for one result object."""

    if isinstance(result, GenePCAResult):
        paths: list[Path] = [plot_gene_pca_variance(result, output_dir)]
        for component_index in range(1, result.variance_table.shape[0] + 1):
            brain_path = plot_brain_volume_map(
                result.regional_scores,
                atlas_id=result.atlas_id,
                value_column=f"PC{component_index}",
                title=f"Brain map for PC{component_index}",
                output_path=output_dir / "plots" / f"gene_pca_pc{component_index}_brain.png",
            )
            if brain_path is not None:
                paths.append(brain_path)
            _append_cortical_plots(
                paths,
                table=result.regional_scores,
                atlas_id=result.atlas_id,
                value_column=f"PC{component_index}",
                title=f"Cortical map for PC{component_index}",
                output_path=output_dir / "plots" / f"gene_pca_pc{component_index}_cortex.png",
            )
            paths.append(plot_gene_pca_regional_component(result, output_dir, component_index))
            paths.append(plot_gene_pca_loadings(result, output_dir, component_index))
        return paths

    if isinstance(result, GEDARResult):
        paths: list[Path] = []
        if result.direction == "split":
            for mode, label in (("up", "Up"), ("down", "Down")):
                score_column = f"score_{mode}_z"
                if score_column not in result.regional_scores.columns:
                    continue
                if np.isnan(result.regional_scores[score_column].to_numpy(dtype=float)).all():
                    continue
                paths.append(
                    plot_gedar_regional_scores(
                        result,
                        output_dir,
                        value_column=score_column,
                        title=f"GEDAR {label.lower()} regional score",
                        filename=f"gedar_{mode}_scores.png",
                    )
                )
                brain_path = plot_brain_volume_map(
                    result.regional_scores,
                    atlas_id=result.atlas_id,
                    value_column=score_column,
                    title=f"GEDAR {label.lower()} brain map",
                    output_path=output_dir / "plots" / f"gedar_{mode}_brain.png",
                )
                if brain_path is not None:
                    paths.append(brain_path)
                _append_cortical_plots(
                    paths,
                    table=result.regional_scores,
                    atlas_id=result.atlas_id,
                    value_column=score_column,
                    title=f"GEDAR {label.lower()} cortical map",
                    output_path=output_dir / "plots" / f"gedar_{mode}_cortex.png",
                )
                weights_path = plot_gedar_weights(
                    result,
                    output_dir,
                    weight_column=f"weight_{mode}",
                    selected_column=f"selected_{mode}",
                    title=f"GEDAR {label.lower()} gene weights",
                    filename=f"gedar_{mode}_weights.png",
                )
                if weights_path is not None:
                    paths.append(weights_path)
        else:
            paths = [plot_gedar_regional_scores(result, output_dir)]
            brain_path = plot_brain_volume_map(
                result.regional_scores,
                atlas_id=result.atlas_id,
                value_column="score_z",
                title="GEDAR brain map",
                output_path=output_dir / "plots" / "gedar_brain.png",
            )
            if brain_path is not None:
                paths.append(brain_path)
            _append_cortical_plots(
                paths,
                table=result.regional_scores,
                atlas_id=result.atlas_id,
                value_column="score_z",
                title="GEDAR cortical map",
                output_path=output_dir / "plots" / "gedar_cortex.png",
            )
            weights_path = plot_gedar_weights(result, output_dir)
            if weights_path is not None:
                paths.append(weights_path)
        return paths

    paths: list[Path] = []
    brain_path = plot_brain_volume_map(
        result.regional_values,
        atlas_id=result.metadata.atlas_id,
        value_column="value",
        title="Regional brain map",
        output_path=output_dir / "plots" / "regional_values_brain.png",
    )
    if brain_path is not None:
        paths.append(brain_path)
    cortical_values = result.regional_values.assign(
        value_z=zscore_for_plot(result.regional_values["value"].to_numpy(dtype=float, copy=False))
    )
    _append_cortical_plots(
        paths,
        table=cortical_values,
        atlas_id=result.metadata.atlas_id,
        value_column="value_z",
        title="Cortical regional map (z-scored)",
        output_path=output_dir / "plots" / "regional_values_cortex.png",
    )
    if isinstance(result, CorrelationResult):
        paths.append(plot_correlation_ranking(result.gene_table, output_dir))
        paths.append(plot_correlation_distribution(result.gene_table, output_dir))
        if result.gsea_table is not None:
            gsea_path = plot_gsea_dotplot(
                result.gsea_table,
                output_dir / "plots" / "gsea_corr_dotplot.png",
                title="Top GSEA terms",
            )
            if gsea_path is not None:
                paths.append(gsea_path)
        if result.ora_tables is not None:
            ora_path = plot_ora_heatmap(
                result.ora_tables,
                output_dir / "plots" / "ora_corr_heatmap.png",
                title="ORA up/down heatmap",
            )
            if ora_path is not None:
                paths.append(ora_path)
        return paths

    paths.extend(plot_pls_variance(result, output_dir))
    for component in result.components:
        paths.append(plot_pls_component(component, output_dir))
        if component.gsea_table is not None:
            gsea_path = plot_gsea_dotplot(
                component.gsea_table,
                output_dir / "plots" / f"gsea_pls{component.index}_dotplot.png",
                title=f"PLS component {component.index} GSEA",
            )
            if gsea_path is not None:
                paths.append(gsea_path)
        if component.ora_tables is not None:
            ora_path = plot_ora_heatmap(
                component.ora_tables,
                output_dir / "plots" / f"ora_pls{component.index}_heatmap.png",
                title=f"PLS component {component.index} ORA up/down heatmap",
            )
            if ora_path is not None:
                paths.append(ora_path)
    return paths
