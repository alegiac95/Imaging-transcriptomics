from __future__ import annotations

from importlib import import_module

__version__ = "2.0.0"

from .api import run_analysis, run_corr, run_gene, run_pls
from .atlas_registry import atlas_table, describe_atlas, get_atlas, list_atlases
from .config import RunConfig, build_run_config
from .exceptions import (
    AtlasAssetError,
    AtlasError,
    ConfigurationError,
    ImagingTranscriptomicsError,
    InputAlignmentError,
    InputDataError,
    NullModelError,
    PlottingUnavailableError,
)
from .gene_expression import load_expression_frame, load_gene_labels, select_atlas_data
from .gene_query import run_gene as run_gene_query
from .gedar import run_gedar
from .gene_pca import run_gene_pca
from .models import (
    AnalysisMetadata,
    AtlasSelection,
    AtlasSpec,
    CorrelationResult,
    ExtractedScan,
    GEDARResult,
    GeneQueryResult,
    GenePCAResult,
    PLSComponentResult,
    PLSResult,
)
from .scan import extract_scan_data, regional_values_frame

__all__ = [
    "__version__",
    "AnalysisMetadata",
    "AtlasAssetError",
    "AtlasError",
    "AtlasSelection",
    "AtlasSpec",
    "ConfigurationError",
    "CorrelationResult",
    "GEDARResult",
    "ExtractedScan",
    "GeneQueryResult",
    "GenePCAResult",
    "ImagingTranscriptomicsError",
    "InputAlignmentError",
    "InputDataError",
    "NullModelError",
    "PlottingUnavailableError",
    "RunConfig",
    "PLSComponentResult",
    "PLSResult",
    "atlas_table",
    "build_run_config",
    "describe_atlas",
    "extract_scan_data",
    "get_atlas",
    "list_atlases",
    "load_expression_frame",
    "load_gene_labels",
    "regional_values_frame",
    "run_analysis",
    "run_corr",
    "run_gene",
    "run_gene_query",
    "run_gedar",
    "run_gene_pca",
    "run_pls",
    "select_atlas_data",
    "build_expression_assets",
]


def __getattr__(name):
    if name == "build_expression_assets":
        return getattr(import_module(".build_atlas", __name__), "build_expression_assets")
    raise AttributeError(name)


def __dir__():
    return sorted(set(globals()) | set(__all__))
