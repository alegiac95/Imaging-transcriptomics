"""Helpers for extracting regional values from volumetric scan inputs."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.processing import resample_from_to

from ..exceptions import AtlasAssetError, InputAlignmentError
from ..models import AtlasSelection


def region_means(image_data: np.ndarray, atlas_data: np.ndarray, atlas_ids: np.ndarray) -> np.ndarray:
    """Compute parcel means from a volumetric image and integer atlas labels."""

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
    """Extract regional means from volumetric data on a packaged atlas grid."""

    atlas_path = selection.atlas.volume_path(resolution)
    if atlas_path is None:
        raise AtlasAssetError(
            f"Atlas '{selection.atlas.id}' does not have a packaged volumetric parcellation for this resolution."
        )
    atlas_data = nib.load(atlas_path).get_fdata()
    atlas_ids = selection.labels["id"].astype(int).to_numpy()
    return region_means(image_data, atlas_data, atlas_ids)


def grid_matches(image: nib.Nifti1Image, atlas_image: nib.Nifti1Image) -> bool:
    """Return True when a NIfTI image already matches an atlas grid."""

    return image.shape[:3] == atlas_image.shape[:3] and np.allclose(
        image.affine,
        atlas_image.affine,
        atol=1e-3,
    )


def atlas_volume_image(selection: AtlasSelection, resolution: str) -> nib.Nifti1Image:
    """Load a packaged volumetric atlas image for a requested resolution."""

    atlas_path = selection.atlas.volume_path(resolution)
    if atlas_path is None:
        raise AtlasAssetError(
            f"Atlas '{selection.atlas.id}' does not have a packaged volumetric parcellation for {resolution}."
        )
    return nib.load(atlas_path)


def matching_atlas_grid(
    image: nib.Nifti1Image,
    selection: AtlasSelection,
) -> tuple[str, nib.Nifti1Image] | None:
    """Return the packaged atlas grid matching the input image, if any."""

    for resolution in ("1mm", "2mm"):
        atlas_path = selection.atlas.volume_path(resolution)
        if atlas_path is None:
            continue
        atlas_image = nib.load(atlas_path)
        if grid_matches(image, atlas_image):
            return resolution, atlas_image
    return None


def preferred_resolution(image: nib.Nifti1Image, selection: AtlasSelection) -> str:
    """Pick the closest packaged atlas resolution for direct resampling."""

    zoom_mean = float(np.mean(np.abs(image.header.get_zooms()[:3])))
    preferred = "1mm" if zoom_mean <= 1.5 else "2mm"
    if selection.atlas.volume_path(preferred) is not None:
        return preferred
    fallback = "2mm" if preferred == "1mm" else "1mm"
    if selection.atlas.volume_path(fallback) is None:
        raise AtlasAssetError(
            f"Atlas '{selection.atlas.id}' does not have a packaged volumetric parcellation for direct extraction."
        )
    return fallback


def native_space_error(image: nib.Nifti1Image, selection: AtlasSelection) -> InputAlignmentError:
    """Create a detailed error for native-space images passed as atlas maps."""

    shape = image.shape[:3]
    zooms = tuple(round(float(val), 3) for val in image.header.get_zooms()[:3])
    return InputAlignmentError(
        "Input NIfTI is not aligned to the packaged MNI atlas grids for "
        f"atlas '{selection.atlas.id}'. Got shape={shape}, zooms={zooms}. "
        "This command does not register native-space subject T1w images into MNI152. "
        "Register the image to MNI152 first, then rerun with `--space MNI152`, "
        "or provide a precomputed regional vector."
    )


def extract_with_atlas_image(
    image_data: np.ndarray,
    atlas_image: nib.Nifti1Image,
    selection: AtlasSelection,
) -> np.ndarray:
    """Extract parcel means using a preloaded atlas image."""

    atlas_ids = selection.labels["id"].astype(int).to_numpy()
    return region_means(image_data, atlas_image.get_fdata(), atlas_ids)


def parcellate_nifti_direct(
    scan_path: Path,
    selection: AtlasSelection,
    *,
    allow_resample: bool,
) -> np.ndarray:
    """Parcellate a volumetric image directly on a packaged MNI atlas grid."""

    image = nib.load(scan_path)
    matched = matching_atlas_grid(image, selection)
    if matched is not None:
        _, atlas_image = matched
        return extract_with_atlas_image(image.get_fdata(), atlas_image, selection)

    if not allow_resample:
        raise native_space_error(image, selection)

    resolution = preferred_resolution(image, selection)
    atlas_image = atlas_volume_image(selection, resolution)
    resampled = resample_from_to(image, atlas_image, order=1)
    return extract_with_atlas_image(resampled.get_fdata(), atlas_image, selection)
