"""Shared helpers for dispatching supported scan input types."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np

from ..exceptions import InputDataError
from ..models import AtlasSelection, ExtractedScan, SourceKind


TEXT_SUFFIXES = {".txt", ".tsv", ".csv"}
NIFTI_SUFFIXES = {".nii", ".nii.gz"}
GIFTI_SUFFIXES = {".gii", ".gii.gz", ".shape.gii", ".func.gii"}
_KNOWN_SUFFIXES = tuple(
    sorted(TEXT_SUFFIXES | NIFTI_SUFFIXES | GIFTI_SUFFIXES, key=len, reverse=True)
)


def recognized_suffix(path: Path) -> str:
    """Return the recognized file suffix from the end of one filename."""

    lower_name = path.name.lower()
    for suffix in _KNOWN_SUFFIXES:
        if lower_name.endswith(suffix):
            return suffix
    return path.suffix.lower()


def make_extracted_scan(
    selection: AtlasSelection,
    values: np.ndarray | Iterable[float],
    *,
    source: str,
    source_space: str | None,
    source_kind: SourceKind,
) -> ExtractedScan:
    """Validate a regional vector and wrap it in the extracted-scan record."""

    vector = np.asarray(values, dtype=float).reshape(-1)
    if vector.shape[0] != selection.n_regions:
        raise InputDataError(
            f"Expected {selection.n_regions} regional values for the selected atlas subset, got {vector.shape[0]}."
        )
    return ExtractedScan(
        values=vector,
        selection=selection,
        source=source,
        source_space=source_space,
        source_kind=source_kind,
    )


def load_tabular_vector(path: Path) -> np.ndarray:
    """Read a one-column text or CSV file as a regional vector."""

    data = np.loadtxt(path, delimiter="," if path.suffix == ".csv" else None)
    if data.ndim > 1:
        if 1 in data.shape:
            data = data.reshape(-1)
        else:
            raise InputDataError(
                f"Expected a one-column vector in {path}, got array with shape {data.shape}."
            )
    return np.asarray(data, dtype=float).reshape(-1)


def surface_inputs(data, input_rh) -> tuple[str, str] | None:
    """Normalize left/right surface inputs into a two-path tuple."""

    if isinstance(data, (list, tuple)) and len(data) == 2:
        return str(Path(data[0])), str(Path(data[1]))
    if input_rh is not None:
        return str(Path(data)), str(Path(input_rh))
    return None
