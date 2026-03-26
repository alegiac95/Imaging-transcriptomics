from __future__ import annotations

from pathlib import Path
from typing import Iterable

import nibabel as nib
import numpy as np
import pandas as pd
from nibabel.processing import resample_from_to

from ._compat import suppress_pkg_resources_deprecation
from .gene_expression import select_atlas_data
from .models import AtlasSelection, ExtractedScan, HemisphereMode, RegionScope
from .surfaces import load_surface_parcellation


_TEXT_SUFFIXES = {".txt", ".tsv", ".csv"}
_NIFTI_SUFFIXES = {".nii", ".nii.gz"}
_GIFTI_SUFFIXES = {".gii", ".gii.gz", ".shape.gii", ".func.gii"}


def _suffix(path: Path) -> str:
    return "".join(path.suffixes) if path.suffix == ".gz" else path.suffix


def _make_extracted_scan(
    selection: AtlasSelection,
    values: np.ndarray | Iterable[float],
    *,
    source: str,
    source_space: str | None,
    source_kind: str,
) -> ExtractedScan:
    vector = np.asarray(values, dtype=float).reshape(-1)
    if vector.shape[0] != selection.n_regions:
        raise ValueError(
            f"Expected {selection.n_regions} regional values for the selected atlas subset, got {vector.shape[0]}."
        )
    return ExtractedScan(
        values=vector,
        selection=selection,
        source=source,
        source_space=source_space,
        source_kind=source_kind,
    )


def _load_tabular_vector(path: Path) -> np.ndarray:
    data = np.loadtxt(path, delimiter="," if path.suffix == ".csv" else None)
    if data.ndim > 1:
        if 1 in data.shape:
            data = data.reshape(-1)
        else:
            raise ValueError(
                f"Expected a one-column vector in {path}, got array with shape {data.shape}."
            )
    return np.asarray(data, dtype=float).reshape(-1)


def _region_means(image_data: np.ndarray, atlas_data: np.ndarray, atlas_ids: np.ndarray) -> np.ndarray:
    atlas_flat = np.asarray(atlas_data, dtype=np.int32).reshape(-1)
    image_flat = np.asarray(image_data, dtype=float).reshape(-1)
    mask = np.isfinite(image_flat) & (atlas_flat > 0)
    if not np.any(mask):
        return np.full(atlas_ids.shape[0], np.nan, dtype=float)

    labels = atlas_flat[mask]
    values = image_flat[mask]
    sums = np.bincount(labels, weights=values)
    counts = np.bincount(labels)
    means = np.full(atlas_ids.shape[0], np.nan, dtype=float)
    valid = atlas_ids < counts.shape[0]
    if np.any(valid):
        valid_ids = atlas_ids[valid]
        valid_counts = counts[valid_ids]
        valid_means = np.full(valid_ids.shape[0], np.nan, dtype=float)
        nonzero = valid_counts > 0
        valid_means[nonzero] = sums[valid_ids[nonzero]] / valid_counts[nonzero]
        means[valid] = valid_means
    return means


def extract_volume_values(
    image_data: np.ndarray,
    selection: AtlasSelection,
    *,
    resolution: str,
) -> np.ndarray:
    atlas_path = selection.atlas.volume_path(resolution)
    if atlas_path is None:
        raise FileNotFoundError(
            f"Atlas '{selection.atlas.id}' does not have a packaged volumetric parcellation for this resolution."
        )
    atlas_data = nib.load(atlas_path).get_fdata()
    atlas_ids = selection.labels["id"].astype(int).to_numpy()
    return _region_means(image_data, atlas_data, atlas_ids)


def _grid_matches(image: nib.Nifti1Image, atlas_image: nib.Nifti1Image) -> bool:
    return image.shape[:3] == atlas_image.shape[:3] and np.allclose(
        image.affine,
        atlas_image.affine,
        atol=1e-3,
    )


def _atlas_volume_image(selection: AtlasSelection, resolution: str) -> nib.Nifti1Image:
    atlas_path = selection.atlas.volume_path(resolution)
    if atlas_path is None:
        raise FileNotFoundError(
            f"Atlas '{selection.atlas.id}' does not have a packaged volumetric parcellation for {resolution}."
        )
    return nib.load(atlas_path)


def _matching_atlas_grid(
    image: nib.Nifti1Image,
    selection: AtlasSelection,
) -> tuple[str, nib.Nifti1Image] | None:
    for resolution in ("1mm", "2mm"):
        atlas_path = selection.atlas.volume_path(resolution)
        if atlas_path is None:
            continue
        atlas_image = nib.load(atlas_path)
        if _grid_matches(image, atlas_image):
            return resolution, atlas_image
    return None


def _preferred_resolution(image: nib.Nifti1Image, selection: AtlasSelection) -> str:
    zoom_mean = float(np.mean(np.abs(image.header.get_zooms()[:3])))
    preferred = "1mm" if zoom_mean <= 1.5 else "2mm"
    if selection.atlas.volume_path(preferred) is not None:
        return preferred
    fallback = "2mm" if preferred == "1mm" else "1mm"
    if selection.atlas.volume_path(fallback) is None:
        raise FileNotFoundError(
            f"Atlas '{selection.atlas.id}' does not have a packaged volumetric parcellation for direct extraction."
        )
    return fallback


def _native_space_error(image: nib.Nifti1Image, selection: AtlasSelection) -> ValueError:
    shape = image.shape[:3]
    zooms = tuple(round(float(val), 3) for val in image.header.get_zooms()[:3])
    return ValueError(
        "Input NIfTI is not aligned to the packaged MNI atlas grids for "
        f"atlas '{selection.atlas.id}'. Got shape={shape}, zooms={zooms}. "
        "This command does not register native-space subject T1w images into MNI152. "
        "Register the image to MNI152 first, then rerun with `--space MNI152`, "
        "or provide a precomputed regional vector."
    )


def _extract_with_atlas_image(
    image_data: np.ndarray,
    atlas_image: nib.Nifti1Image,
    selection: AtlasSelection,
) -> np.ndarray:
    atlas_ids = selection.labels["id"].astype(int).to_numpy()
    return _region_means(image_data, atlas_image.get_fdata(), atlas_ids)


def _parcellate_nifti_direct(
    scan_path: Path,
    selection: AtlasSelection,
    *,
    allow_resample: bool,
) -> np.ndarray:
    image = nib.load(scan_path)
    matched = _matching_atlas_grid(image, selection)
    if matched is not None:
        _, atlas_image = matched
        return _extract_with_atlas_image(image.get_fdata(), atlas_image, selection)

    if not allow_resample:
        raise _native_space_error(image, selection)

    resolution = _preferred_resolution(image, selection)
    atlas_image = _atlas_volume_image(selection, resolution)
    resampled = resample_from_to(image, atlas_image, order=1)
    return _extract_with_atlas_image(resampled.get_fdata(), atlas_image, selection)


def _import_neuromaps():
    try:
        with suppress_pkg_resources_deprecation():
            from neuromaps.parcellate import Parcellater
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "neuromaps is required for cross-space resampling or surface inputs. "
            "Install imaging-transcriptomics[maps]."
        ) from exc
    return Parcellater


def _parcellate_with_neuromaps(
    data,
    selection: AtlasSelection,
    *,
    source_space: str,
    hemi: str | None = None,
) -> np.ndarray:
    Parcellater = _import_neuromaps()
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
        raise ValueError(
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
        raise ValueError(
            f"Parcellated data returned {values.shape[0]} values, which is too short for atlas ids up to {int(atlas_ids.max())}."
        )
    return values[atlas_ids - 1]


def _surface_inputs(data, input_rh) -> tuple[str, str] | None:
    if isinstance(data, (list, tuple)) and len(data) == 2:
        return str(Path(data[0])), str(Path(data[1]))
    if input_rh is not None:
        return str(Path(data)), str(Path(input_rh))
    return None


def extract_scan_data(
    data,
    atlas: str = "dk",
    hemisphere: HemisphereMode = "left",
    regions: RegionScope = "all",
    source_space: str | None = None,
    input_rh=None,
    prefer_neuromaps: bool = True,
) -> ExtractedScan:
    selection = select_atlas_data(atlas=atlas, hemisphere=hemisphere, regions=regions)

    if isinstance(data, np.ndarray):
        return _make_extracted_scan(
            selection,
            data,
            source="array",
            source_space=source_space,
            source_kind="vector",
        )

    surface_inputs = _surface_inputs(data, input_rh)
    if surface_inputs is not None:
        values = _parcellate_with_neuromaps(
            surface_inputs,
            selection,
            source_space=source_space or "fsaverage",
        )
        return _make_extracted_scan(
            selection,
            values,
            source=",".join(surface_inputs),
            source_space=source_space or "fsaverage",
            source_kind="surface",
        )

    path = Path(data)
    suffix = _suffix(path)
    if suffix in _TEXT_SUFFIXES:
        return _make_extracted_scan(
            selection,
            _load_tabular_vector(path),
            source=str(path),
            source_space=source_space,
            source_kind="vector",
        )
    if suffix in _NIFTI_SUFFIXES:
        if source_space is None:
            values = _parcellate_nifti_direct(path, selection, allow_resample=False)
            resolved_space = "MNI152"
        elif source_space == "MNI152":
            values = _parcellate_nifti_direct(path, selection, allow_resample=True)
            resolved_space = "MNI152"
        elif prefer_neuromaps:
            values = _parcellate_with_neuromaps(path.as_posix(), selection, source_space=source_space)
            resolved_space = source_space
        else:
            raise ValueError(
                "Non-MNI volumetric data require neuromaps-based resampling. Pass source_space and install neuromaps."
            )
        return _make_extracted_scan(
            selection,
            values,
            source=str(path),
            source_space=resolved_space,
            source_kind="volume",
        )
    if suffix in _GIFTI_SUFFIXES:
        raise ValueError("Surface inputs require both left and right hemisphere files.")
    raise ValueError(
        f"Unsupported input type for '{data}'. Expected a vector, NIfTI, text table, or left/right GIFTI pair."
    )


def regional_values_frame(extracted: ExtractedScan) -> pd.DataFrame:
    return extracted.labels.assign(value=extracted.values)
