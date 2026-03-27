from __future__ import annotations


class ImagingTranscriptomicsError(Exception):
    """Base exception for package-specific failures."""


class ConfigurationError(ValueError, ImagingTranscriptomicsError):
    """Raised when user-facing configuration values are invalid."""


class AtlasError(ImagingTranscriptomicsError):
    """Base exception for atlas lookup and atlas-asset failures."""


class AtlasAssetError(FileNotFoundError, AtlasError):
    """Raised when required atlas assets are missing."""


class InputDataError(ValueError, ImagingTranscriptomicsError):
    """Raised when an input file or vector cannot be interpreted."""


class InputAlignmentError(InputDataError):
    """Raised when an image exists but is not aligned to the expected space."""


class NullModelError(RuntimeError, ImagingTranscriptomicsError):
    """Raised when null-model generation fails."""


class PlottingUnavailableError(ImportError, ImagingTranscriptomicsError):
    """Raised when optional plotting dependencies are unavailable."""
