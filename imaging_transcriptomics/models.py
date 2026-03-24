from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd


HemisphereMode = Literal["left", "both"]
RegionScope = Literal["all", "cort", "cort+sub"]
AnalysisMethod = Literal["corr", "pls"]
NullMethod = Literal["auto", "vasa", "alexander_bloch", "moran", "random"]


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
    volume_1mm_path: Path | None = None
    volume_2mm_path: Path | None = None
    lh_annot_path: Path | None = None
    rh_annot_path: Path | None = None
    notes: str | None = None

    def volume_path(self, resolution: str) -> Path | None:
        if resolution == "1mm":
            return self.volume_1mm_path
        if resolution == "2mm":
            return self.volume_2mm_path
        raise ValueError(f"Unsupported atlas resolution: {resolution}")


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
    source_kind: str = "vector"

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
    source_kind: str
    source_space: str | None
    n_permutations: int
    null_method: str = "auto"
    geneset: str | None = None
    n_components: int | None = None


@dataclass(frozen=True)
class CorrelationResult:
    """Returned by :func:`run_corr`."""

    metadata: AnalysisMetadata
    regional_values: pd.DataFrame
    gene_table: pd.DataFrame
    gsea_table: pd.DataFrame | None = None
    output_dir: Path | None = None


@dataclass(frozen=True)
class PLSComponentResult:
    """One component from a PLS run."""

    index: int
    explained_variance: float
    p_value: float
    gene_table: pd.DataFrame
    gsea_table: pd.DataFrame | None = None


@dataclass(frozen=True)
class PLSResult:
    """Returned by :func:`run_pls`."""

    metadata: AnalysisMetadata
    regional_values: pd.DataFrame
    components: tuple[PLSComponentResult, ...]
    cumulative_variance: np.ndarray
    output_dir: Path | None = None
