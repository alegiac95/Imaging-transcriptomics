from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd


HemisphereMode = Literal["left", "both"]
RegionScope = Literal["default", "all", "cort", "cort+sub"]
AnalysisMethod = Literal["corr", "pls", "gene-pca", "gedar"]
NullMethod = Literal["auto", "vasa", "alexander_bloch", "moran", "random"]
SourceKind = Literal["vector", "surface", "volume"]
RankMode = Literal["ascending", "descending"]
GEDARDirection = Literal["combined", "up", "down", "split"]
ExpressionNormalization = Literal["zscore", "none"]
WeightNormalization = Literal["none", "zscore", "unit"]


@dataclass(frozen=True)
class AtlasSpec:
    """Atlas metadata used by the v2 functional API."""

    id: str
    label: str
    description: str
    family: str
    volumetric_space: str
    supported_spaces: tuple[str, ...]
    packaged: bool
    buildable: bool
    default_hemisphere: HemisphereMode
    has_subcortex: bool
    n_regions_left: int
    n_regions_both: int
    labels_path: Path | None = None
    expression_path: Path | None = None
    gene_labels_path: Path | None = None
    volume_1mm_path: Path | None = None
    volume_2mm_path: Path | None = None
    lh_annot_path: Path | None = None
    rh_annot_path: Path | None = None
    lh_surface_path: Path | None = None
    rh_surface_path: Path | None = None
    surface_space: str | None = None
    surface_density: str | None = None
    geometry_lh_path: Path | None = None
    geometry_rh_path: Path | None = None
    notes: str | None = None

    def volume_path(self, resolution: str) -> Path | None:
        if resolution == "1mm":
            return self.volume_1mm_path
        if resolution == "2mm":
            return self.volume_2mm_path
        raise ValueError(f"Unsupported atlas resolution: {resolution}")

    @property
    def surface_paths(self) -> tuple[Path, Path] | None:
        if self.lh_surface_path is not None and self.rh_surface_path is not None:
            return self.lh_surface_path, self.rh_surface_path
        if self.lh_annot_path is not None and self.rh_annot_path is not None:
            return self.lh_annot_path, self.rh_annot_path
        return None

    @property
    def surface_geometry(self) -> tuple[Path, Path] | None:
        if self.geometry_lh_path is not None and self.geometry_rh_path is not None:
            return self.geometry_lh_path, self.geometry_rh_path
        return None


@dataclass(frozen=True)
class AtlasSelection:
    """Concrete atlas subset used for one analysis run."""

    atlas: AtlasSpec
    hemisphere: HemisphereMode
    regions: RegionScope
    labels: pd.DataFrame
    expression: pd.DataFrame
    gene_labels: np.ndarray

    @property
    def n_regions(self) -> int:
        return int(self.labels.shape[0])

    @property
    def region_names(self) -> list[str]:
        return self.labels["label"].astype(str).tolist()


@dataclass(frozen=True)
class ExtractedScan:
    """Parcellated scan values aligned to one atlas selection."""

    values: np.ndarray
    selection: AtlasSelection
    source: str
    source_space: str | None = None
    source_kind: SourceKind = "vector"

    @property
    def labels(self) -> pd.DataFrame:
        return self.selection.labels

    @property
    def region_names(self) -> list[str]:
        return self.selection.region_names


@dataclass(frozen=True)
class AnalysisMetadata:
    """Small bundle of run settings that should travel with the results."""

    method: AnalysisMethod
    atlas_id: str
    atlas_label: str
    hemisphere: HemisphereMode
    regions: RegionScope
    source: str
    source_kind: SourceKind
    source_space: str | None
    n_permutations: int
    null_method: NullMethod = "auto"
    geneset: str | None = None
    geneset_organism: str | None = None
    ora_p_threshold: float | None = None
    n_components: int | None = None


@dataclass(frozen=True)
class CorrelationResult:
    """Returned by :func:`run_corr`."""

    metadata: AnalysisMetadata
    regional_values: pd.DataFrame
    gene_table: pd.DataFrame
    gsea_table: pd.DataFrame | None = None
    ora_tables: dict[str, pd.DataFrame] | None = None
    output_dir: Path | None = None


@dataclass(frozen=True)
class PLSComponentResult:
    """One component from a PLS run."""

    index: int
    explained_variance: float
    p_value: float
    gene_table: pd.DataFrame
    gsea_table: pd.DataFrame | None = None
    ora_tables: dict[str, pd.DataFrame] | None = None


@dataclass(frozen=True)
class PLSResult:
    """Returned by :func:`run_pls`."""

    metadata: AnalysisMetadata
    regional_values: pd.DataFrame
    components: tuple[PLSComponentResult, ...]
    cumulative_variance: np.ndarray
    output_dir: Path | None = None


@dataclass(frozen=True)
class GenePCAResult:
    """Returned by :func:`run_gene_pca`."""

    atlas_id: str
    atlas_label: str
    hemisphere: HemisphereMode
    regions: RegionScope
    requested_genes: tuple[str, ...]
    regional_scores: pd.DataFrame
    gene_loadings: pd.DataFrame
    variance_table: pd.DataFrame
    matched_genes: tuple[str, ...]
    brain_filtered_genes: tuple[str, ...]
    missing_genes: tuple[str, ...]
    output_dir: Path | None = None


@dataclass(frozen=True)
class GEDARResult:
    """Returned by :func:`run_gedar`."""

    atlas_id: str
    atlas_label: str
    hemisphere: HemisphereMode
    regions: RegionScope
    requested_genes: tuple[str, ...]
    regional_scores: pd.DataFrame
    gene_table: pd.DataFrame
    excluded_table: pd.DataFrame
    matched_genes: tuple[str, ...]
    missing_genes: tuple[str, ...]
    weights_source: str
    gene_column: str
    weight_column: str
    rank_column: str | None
    rank_mode: RankMode
    direction: GEDARDirection
    normalize_expression: ExpressionNormalization
    normalize_weights: WeightNormalization
    top_percent: float | None
    top_n: int | None
    p_threshold: float | None
    output_dir: Path | None = None
