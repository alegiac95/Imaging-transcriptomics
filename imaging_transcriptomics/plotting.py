"""Public plotting facade for the imaging transcriptomics toolbox.

The plotting implementation lives in smaller modules under
``imaging_transcriptomics.outputs``. This facade keeps the public import path
stable for the rest of the package and for downstream users.
"""

from __future__ import annotations

from .outputs.brain import (
    load_surface_parcellation as _load_surface_parcellation,
    plot_brain_volume_map,
    plot_cortical_surface_map,
    plot_cortical_surface_map_brainspace,
    surface_value_frames as _surface_value_frames,
    vertex_values_for_hemisphere as _vertex_values_for_hemisphere,
)
from .outputs.bundle import save_result_plots
from .outputs.common import (
    matplotlib_backend as _matplotlib,
    regional_colors as _regional_colors,
    safe_neglog10 as _safe_neglog10,
    shorten_labels as _shorten_labels,
    style_axes as _style_axes,
    zscore_for_plot as _zscore_for_plot,
)
from .outputs.enrichment import (
    ensemble_dot_frame as _ensemble_dot_frame,
    ora_heatmap_frame as _ora_heatmap_frame,
    plot_ensemble_dotplot,
    plot_gsea_dotplot,
    plot_ora_heatmap,
)
from .outputs.gene_tables import (
    plot_correlation_distribution,
    plot_correlation_ranking,
    plot_gene_query_distribution,
    plot_gene_query_matrix,
    plot_gene_query_ranking,
    plot_gene_pca_loadings,
    plot_gene_pca_regional_component,
    plot_gene_pca_variance,
    plot_gedar_regional_scores,
    plot_gedar_weights,
    plot_pls_component,
    plot_pls_variance,
    plot_region_profile,
    plot_regional_values,
)

__all__ = [
    "_load_surface_parcellation",
    "_ensemble_dot_frame",
    "_matplotlib",
    "_ora_heatmap_frame",
    "_regional_colors",
    "_safe_neglog10",
    "_shorten_labels",
    "_style_axes",
    "_surface_value_frames",
    "_vertex_values_for_hemisphere",
    "_zscore_for_plot",
    "plot_brain_volume_map",
    "plot_correlation_distribution",
    "plot_correlation_ranking",
    "plot_gene_query_distribution",
    "plot_gene_query_matrix",
    "plot_gene_query_ranking",
    "plot_cortical_surface_map",
    "plot_cortical_surface_map_brainspace",
    "plot_gene_pca_loadings",
    "plot_gene_pca_regional_component",
    "plot_gene_pca_variance",
    "plot_gedar_regional_scores",
    "plot_gedar_weights",
    "plot_ensemble_dotplot",
    "plot_gsea_dotplot",
    "plot_ora_heatmap",
    "plot_pls_component",
    "plot_pls_variance",
    "plot_region_profile",
    "plot_regional_values",
    "save_result_plots",
]
