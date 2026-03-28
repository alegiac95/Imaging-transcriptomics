from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from .._compat import suppress_pkg_resources_deprecation
from ..atlas_registry import get_atlas
from .common import BACKGROUND, matplotlib_backend, save_figure


def atlas_volume_path(atlas_id: str) -> Path | None:
    """Return the preferred packaged atlas volume used for brain slice plots."""

    atlas = get_atlas(atlas_id)
    return atlas.volume_2mm_path or atlas.volume_1mm_path


@lru_cache(maxsize=None)
def fetch_standard_surface_meshes(space: str, density: str) -> tuple[tuple[str, str], tuple[str, str]] | None:
    """Fetch reference surface meshes for atlases without packaged geometry."""

    try:
        with suppress_pkg_resources_deprecation():
            from neuromaps.datasets import fetch_atlas
    except ImportError:
        return None

    atlas = fetch_atlas(space, density, verbose=0)
    inflated = atlas.inflated
    sulc = getattr(atlas, "sulc", None)
    inflated_paths = (str(inflated[0]), str(inflated[1]))
    if sulc is None:
        return inflated_paths, ("", "")
    return inflated_paths, (str(sulc[0]), str(sulc[1]))


def surface_mesh_paths(atlas_id: str) -> tuple[tuple[str, str], tuple[str, str]] | None:
    """Return surface mesh paths for one atlas if cortical plotting is possible."""

    atlas = get_atlas(atlas_id)
    if atlas.surface_geometry is not None:
        inflated = (str(atlas.surface_geometry[0]), str(atlas.surface_geometry[1]))
        return inflated, ("", "")
    if atlas.surface_space is None or atlas.surface_density is None:
        return None
    return fetch_standard_surface_meshes(atlas.surface_space, atlas.surface_density)


def load_surface_mesh(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Load one surface mesh as point and triangle arrays."""

    surface = nib.load(path)
    points = np.asarray(surface.agg_data("pointset"), dtype=float)
    triangles = np.asarray(surface.agg_data("triangle"), dtype=np.int32)
    return points, triangles


def normalize_surface_region_name(name: str) -> str:
    """Normalize atlas label names across FreeSurfer and GIFTI conventions."""

    text = str(name).strip()
    lower = text.lower()
    for prefix in ("ctx-lh-", "ctx-rh-"):
        if lower.startswith(prefix):
            return text[len(prefix) :]
    return text


def load_surface_parcellation(path: str) -> tuple[np.ndarray, dict[int, str]]:
    """Load surface parcel codes together with their normalized region names."""

    surface_path = Path(path)
    if surface_path.suffix == ".annot":
        from nibabel.freesurfer import read_annot

        labels, _, names = read_annot(str(surface_path))
        code_to_name = {
            int(code): normalize_surface_region_name(
                name.decode() if isinstance(name, bytes) else str(name)
            )
            for code, name in enumerate(names)
        }
        return np.asarray(labels, dtype=np.int32), code_to_name

    image = nib.load(str(surface_path))
    data = np.asarray(image.agg_data(), dtype=np.int32)
    code_to_name: dict[int, str] = {}
    label_table = getattr(image, "labeltable", None)
    if label_table is not None:
        for label in label_table.labels:
            code_to_name[int(label.key)] = normalize_surface_region_name(label.label)
    return data, code_to_name


def surface_value_frames(table: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Split one regional table into left/right cortical subsets for surface plots."""

    if not {"hemisphere", "structure"}.issubset(table.columns):
        return {}
    cortical = table.loc[table["structure"].astype(str) == "cortex"].copy()
    if cortical.empty:
        return {}
    frames: dict[str, pd.DataFrame] = {}
    for hemi, label in (("L", "left"), ("R", "right")):
        frame = cortical.loc[cortical["hemisphere"].astype(str) == hemi].copy()
        if not frame.empty:
            frames[label] = frame
    return frames


def vertex_values_for_hemisphere(
    frame: pd.DataFrame,
    *,
    value_column: str,
    label_array: np.ndarray,
    code_to_name: dict[int, str] | None = None,
) -> np.ndarray:
    """Expand regional values to one hemisphere's surface vertices."""

    values = np.full(label_array.shape, np.nan, dtype=float)
    if "label" in frame.columns and code_to_name:
        region_lookup = {
            normalize_surface_region_name(label): value
            for label, value in zip(
                frame["label"].astype(str),
                frame[value_column].to_numpy(dtype=float, copy=False),
                strict=False,
            )
        }
        assigned = 0
        for code, name in code_to_name.items():
            region_value = region_lookup.get(name)
            if region_value is None:
                continue
            values[label_array == code] = region_value
            assigned += 1
        if assigned:
            return values

    region_ids = frame["id"].to_numpy(dtype=int, copy=False)
    region_values = frame[value_column].to_numpy(dtype=float, copy=False)
    for region_id, region_value in zip(region_ids, region_values, strict=False):
        values[label_array == region_id] = region_value
    return values


def surface_face_colors(face_values: np.ndarray, *, cmap, norm) -> np.ndarray:
    """Return RGBA face colors for one array of triangle values."""

    colors = np.asarray(cmap(norm(np.nan_to_num(face_values, nan=0.0))), dtype=float)
    colors[~np.isfinite(face_values)] = np.array([0.90, 0.93, 0.97, 1.0], dtype=float)
    return colors


def triangle_face_values(vertex_values: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Average vertex values onto each triangle face, keeping missing faces blank."""

    tri_values = vertex_values[triangles]
    valid = np.isfinite(tri_values)
    counts = valid.sum(axis=1)
    sums = np.where(valid, tri_values, 0.0).sum(axis=1)
    return np.divide(
        sums,
        counts,
        out=np.full(triangles.shape[0], np.nan, dtype=float),
        where=counts > 0,
    )


def surface_view(ax, coords: np.ndarray, triangles: np.ndarray, vertex_values: np.ndarray, *, title: str, azim: float, cmap, norm) -> None:
    """Render one cortical surface panel."""

    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    face_values = triangle_face_values(vertex_values, triangles)
    face_colors = surface_face_colors(face_values, cmap=cmap, norm=norm)
    surf = ax.plot_trisurf(
        coords[:, 0],
        coords[:, 1],
        coords[:, 2],
        triangles=triangles,
        linewidth=0.0,
        antialiased=False,
        shade=False,
    )
    surf.set_facecolors(face_colors)
    surf.set_edgecolor("none")
    ax.view_init(elev=0, azim=azim)
    ax.set_title(title, fontsize=10)
    ax.set_axis_off()
    extent = np.ptp(coords, axis=0)
    ax.set_box_aspect(tuple(np.maximum(extent, 1e-3)))
    ax.set_facecolor("white")


def plot_cortical_surface_map(
    table: pd.DataFrame,
    *,
    atlas_id: str,
    value_column: str,
    title: str,
    output_path: Path,
) -> Path | None:
    """Render one regional table onto the atlas cortical surface."""

    _, plt = matplotlib_backend()
    atlas = get_atlas(atlas_id)
    if atlas.surface_paths is None:
        return None

    hemi_frames = surface_value_frames(table)
    if not hemi_frames:
        return None

    mesh_paths = surface_mesh_paths(atlas_id)
    if mesh_paths is None:
        return None
    inflated_paths, _ = mesh_paths

    panels: list[tuple[str, str, float]] = []
    if "left" in hemi_frames:
        panels.extend([("left", "L lateral", 180), ("left", "L medial", 0)])
    if "right" in hemi_frames:
        panels.extend([("right", "R medial", 180), ("right", "R lateral", 0)])
    if not panels:
        return None

    finite_values = table[value_column].to_numpy(dtype=float, copy=False)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size == 0:
        return None
    vmax = float(np.nanpercentile(np.abs(finite_values), 98))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = float(np.nanmax(np.abs(finite_values))) if finite_values.size else 1.0
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0

    from matplotlib.colors import TwoSlopeNorm

    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    cmap = plt.get_cmap("RdBu_r")
    fig = plt.figure(figsize=(4.0 * len(panels), 3.4))
    fig.patch.set_facecolor("white")

    hemi_index = {"left": 0, "right": 1}
    for panel_index, (hemi, panel_title, azim) in enumerate(panels, start=1):
        ax = fig.add_subplot(1, len(panels), panel_index, projection="3d")
        coords, triangles = load_surface_mesh(inflated_paths[hemi_index[hemi]])
        label_array, code_to_name = load_surface_parcellation(
            str(atlas.surface_paths[hemi_index[hemi]])
        )
        vertex_values = vertex_values_for_hemisphere(
            hemi_frames[hemi],
            value_column=value_column,
            label_array=label_array,
            code_to_name=code_to_name,
        )
        surface_view(ax, coords, triangles, vertex_values, title=panel_title, azim=azim, cmap=cmap, norm=norm)

    fig.suptitle(title, x=0.04, y=0.98, ha="left", fontweight="bold", fontsize=12)
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])
    colorbar = fig.colorbar(mappable, ax=fig.axes, pad=0.02, shrink=0.72)
    colorbar.set_label(value_column)
    return save_figure(fig, output_path)


def brain_volume_from_table(table: pd.DataFrame, atlas_id: str, value_column: str) -> np.ndarray | None:
    """Project one regional table into a volumetric atlas array."""

    if "id" not in table.columns or value_column not in table.columns:
        return None
    atlas_path = atlas_volume_path(atlas_id)
    if atlas_path is None:
        return None
    atlas_data = np.asarray(nib.load(atlas_path).get_fdata(), dtype=np.int32)
    values = np.full(atlas_data.shape, np.nan, dtype=float)
    ids = table["id"].to_numpy(dtype=int, copy=False)
    score_values = table[value_column].to_numpy(dtype=float, copy=False)
    for region_id, score in zip(ids, score_values, strict=False):
        values[atlas_data == region_id] = score
    if not np.any(np.isfinite(values)):
        return None
    return values


def slice_index(volume: np.ndarray, axis: int) -> int:
    """Choose the most informative slice index for one axis."""

    mask = np.isfinite(volume)
    if not np.any(mask):
        return volume.shape[axis] // 2
    other_axes = tuple(idx for idx in range(3) if idx != axis)
    counts = np.sum(mask, axis=other_axes)
    return int(np.argmax(counts))


def slice_view(volume: np.ndarray, axis: int, index: int) -> np.ndarray:
    """Extract and rotate one orthogonal slice for display."""

    if axis == 0:
        view = volume[index, :, :]
    elif axis == 1:
        view = volume[:, index, :]
    else:
        view = volume[:, :, index]
    return np.rot90(view)


def plot_brain_volume_map(
    table: pd.DataFrame,
    *,
    atlas_id: str,
    value_column: str,
    title: str,
    output_path: Path,
) -> Path | None:
    """Render three orthogonal atlas slice views for one regional table."""

    _, plt = matplotlib_backend()
    volume = brain_volume_from_table(table, atlas_id, value_column)
    if volume is None:
        return None

    finite_values = volume[np.isfinite(volume)]
    if finite_values.size == 0:
        return None
    vmax = float(np.nanpercentile(np.abs(finite_values), 98))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = float(np.nanmax(np.abs(finite_values))) if finite_values.size else 1.0
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0

    from matplotlib.colors import TwoSlopeNorm

    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.2))
    fig.patch.set_facecolor("white")
    slice_specs = [
        ("Sagittal", 0, slice_index(volume, 0)),
        ("Coronal", 1, slice_index(volume, 1)),
        ("Axial", 2, slice_index(volume, 2)),
    ]
    image = None
    for ax, (label, axis, index) in zip(axes, slice_specs, strict=False):
        data = np.ma.masked_invalid(slice_view(volume, axis, index))
        image = ax.imshow(data, cmap="RdBu_r", norm=norm, interpolation="nearest")
        ax.set_title(label, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_facecolor(BACKGROUND)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle(title, x=0.06, y=0.98, ha="left", fontweight="bold", fontsize=12)
    if image is not None:
        colorbar = fig.colorbar(image, ax=axes, pad=0.02, shrink=0.85)
        colorbar.set_label(value_column)
    return save_figure(fig, output_path)
