"""Public scan-loading facade for atlas-aligned regional value extraction."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .exceptions import InputDataError
from .gene_expression import select_atlas_data
from .inputs import (
    GIFTI_SUFFIXES,
    NIFTI_SUFFIXES,
    TEXT_SUFFIXES,
    load_tabular_vector,
    make_extracted_scan,
    parcellate_nifti_direct,
    parcellate_with_neuromaps,
    recognized_suffix,
    surface_inputs,
)
from .models import ExtractedScan, HemisphereMode, RegionScope


def extract_scan_data(
    data,
    atlas: str = "dk",
    hemisphere: HemisphereMode = "left",
    regions: RegionScope = "default",
    source_space: str | None = None,
    input_rh=None,
    prefer_neuromaps: bool = True,
) -> ExtractedScan:
    """Extract atlas-aligned regional values from a supported imaging input.

    Parameters
    ----------
    data
        Input imaging data. Accepted forms are a NumPy vector, a text file
        containing one regional value per row, a NIfTI file, or a left/right
        surface pair.
    atlas
        Atlas identifier registered in the packaged atlas registry.
    hemisphere
        Hemisphere subset to select from the atlas expression data.
    regions
        Region subset to select from the atlas expression data.
    source_space
        Declared standard space of the input image. When omitted for NIfTI
        inputs, the code expects the image to already match a packaged MNI grid.
    input_rh
        Optional right-hemisphere surface file when ``data`` points to the left
        hemisphere file.
    prefer_neuromaps
        Whether to use ``neuromaps`` for supported cross-space resampling.

    Returns
    -------
    ExtractedScan
        Regional values together with the atlas selection and basic source
        metadata.
    """

    selection = select_atlas_data(atlas=atlas, hemisphere=hemisphere, regions=regions)

    if isinstance(data, np.ndarray):
        return make_extracted_scan(
            selection,
            data,
            source="array",
            source_space=source_space,
            source_kind="vector",
        )

    surfaces = surface_inputs(data, input_rh)
    if surfaces is not None:
        values = parcellate_with_neuromaps(
            surfaces,
            selection,
            source_space=source_space or "fsaverage",
        )
        return make_extracted_scan(
            selection,
            values,
            source=",".join(surfaces),
            source_space=source_space or "fsaverage",
            source_kind="surface",
        )

    path = Path(data)
    suffix = recognized_suffix(path)
    if suffix in TEXT_SUFFIXES:
        return make_extracted_scan(
            selection,
            load_tabular_vector(path),
            source=str(path),
            source_space=source_space,
            source_kind="vector",
        )
    if suffix in NIFTI_SUFFIXES:
        if source_space is None:
            values = parcellate_nifti_direct(path, selection, allow_resample=False)
            resolved_space = "MNI152"
        elif source_space == "MNI152":
            values = parcellate_nifti_direct(path, selection, allow_resample=True)
            resolved_space = "MNI152"
        elif prefer_neuromaps:
            values = parcellate_with_neuromaps(path.as_posix(), selection, source_space=source_space)
            resolved_space = source_space
        else:
            raise InputDataError(
                "Non-MNI volumetric data require neuromaps-based resampling. Pass source_space and install neuromaps."
            )
        return make_extracted_scan(
            selection,
            values,
            source=str(path),
            source_space=resolved_space,
            source_kind="volume",
        )
    if suffix in GIFTI_SUFFIXES:
        raise InputDataError("Surface inputs require both left and right hemisphere files.")
    raise InputDataError(
        f"Unsupported input type for '{data}'. Expected a vector, NIfTI, text table, or left/right GIFTI pair."
    )


def regional_values_frame(extracted: ExtractedScan) -> pd.DataFrame:
    """Return regional values merged with atlas label metadata."""

    return extracted.labels.assign(value=extracted.values)
