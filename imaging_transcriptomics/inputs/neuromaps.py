"""Lazy neuromaps-backed scan parcellation helpers."""

from __future__ import annotations

import numpy as np

from .._compat import suppress_pkg_resources_deprecation
from ..exceptions import AtlasAssetError, InputDataError
from ..models import AtlasSelection
from ..surfaces import load_surface_parcellation


def import_neuromaps():
    """Import neuromaps lazily so it stays an optional dependency."""

    try:
        with suppress_pkg_resources_deprecation():
            from neuromaps.parcellate import Parcellater
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "neuromaps is required for cross-space resampling or surface inputs "
            "and is part of the standard imaging-transcriptomics install. "
            "Reinstall the package if it is missing from the current environment."
        ) from exc
    return Parcellater


def parcellate_with_neuromaps(
    data,
    selection: AtlasSelection,
    *,
    source_space: str,
    hemi: str | None = None,
) -> np.ndarray:
    """Parcellate data after neuromaps resampling between standard spaces."""

    Parcellater = import_neuromaps()
    atlas = selection.atlas
    if atlas.volume_1mm_path is not None and source_space == "MNI152" and not isinstance(data, tuple):
        parcellater = Parcellater(
            atlas.volume_1mm_path,
            space="MNI152",
            resampling_target="parcellation",
        )
    elif atlas.surface_paths is not None and atlas.surface_space is not None:
        parcellation = load_surface_parcellation(atlas, "both")
        parcellater = Parcellater(
            parcellation,
            space=atlas.surface_space,
            resampling_target="parcellation",
        )
    else:  # pragma: no cover - depends on optional atlas assets
        raise AtlasAssetError(
            f"Atlas '{atlas.id}' does not have the assets required for neuromaps parcellation."
        )

    values = np.asarray(
        parcellater.fit_transform(
            data,
            source_space,
            ignore_background_data=True,
            hemi=hemi,
        )
    ).reshape(-1)
    atlas_ids = selection.labels["id"].astype(int).to_numpy()
    if values.shape[0] < int(atlas_ids.max()):
        raise InputDataError(
            f"Parcellated data returned {values.shape[0]} values, which is too short for atlas ids up to {int(atlas_ids.max())}."
        )
    return values[atlas_ids - 1]
