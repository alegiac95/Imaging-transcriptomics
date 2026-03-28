"""Internal helpers for parsing and parcellating supported scan inputs."""

from .common import (
    GIFTI_SUFFIXES,
    NIFTI_SUFFIXES,
    TEXT_SUFFIXES,
    load_tabular_vector,
    make_extracted_scan,
    recognized_suffix,
    surface_inputs,
)
from .neuromaps import parcellate_with_neuromaps
from .volume import parcellate_nifti_direct

__all__ = [
    "GIFTI_SUFFIXES",
    "NIFTI_SUFFIXES",
    "TEXT_SUFFIXES",
    "load_tabular_vector",
    "make_extracted_scan",
    "parcellate_nifti_direct",
    "parcellate_with_neuromaps",
    "recognized_suffix",
    "surface_inputs",
]
