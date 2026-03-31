from __future__ import annotations

import importlib.util
import logging
import os
import subprocess
import sys
import tempfile
from functools import lru_cache
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from .._compat import suppress_pkg_resources_deprecation
from ..atlas_registry import get_atlas
from .common import BACKGROUND, MISSING, matplotlib_backend, save_figure

logger = logging.getLogger(__name__)


def atlas_volume_path(atlas_id: str) -> Path | None:
    """Return the preferred packaged atlas volume used for brain slice plots."""

    atlas = get_atlas(atlas_id)
    return atlas.volume_2mm_path or atlas.volume_1mm_path


@lru_cache(maxsize=None)
def fetch_standard_surface_meshes(
    space: str,
    density: str,
    mesh_kind: str = "pial",
) -> tuple[str, str] | None:
    """Fetch reference surface meshes for atlases without packaged geometry.

    Parameters
    ----------
    space, density
        Standard surface space specifiers understood by :mod:`neuromaps`.
    mesh_kind
        Preferred mesh family. ``"pial"`` gives a more anatomical cortical
        surface, while ``"inflated"`` yields a smoother publication-style
        rendering with less visible sulcal folding.
    """

    try:
        with suppress_pkg_resources_deprecation():
            from neuromaps.datasets import fetch_atlas
    except ImportError:
        return None

    atlas = fetch_atlas(space, density, verbose=0)
    candidates: list[object] = []
    if mesh_kind == "inflated":
        candidates.extend(
            [
                getattr(atlas, "inflated", None),
                getattr(atlas, "pial", None),
                getattr(atlas, "white", None),
            ]
        )
    else:
        candidates.extend(
            [
                getattr(atlas, "pial", None),
                getattr(atlas, "inflated", None),
                getattr(atlas, "white", None),
            ]
        )
    mesh = next((candidate for candidate in candidates if candidate is not None), None)
    if mesh is None:
        return None
    return str(mesh[0]), str(mesh[1])


def surface_mesh_paths(atlas_id: str, *, mesh_kind: str = "pial") -> tuple[str, str] | None:
    """Return surface mesh paths for one atlas if cortical plotting is possible."""

    atlas = get_atlas(atlas_id)
    if atlas.surface_geometry is not None:
        return str(atlas.surface_geometry[0]), str(atlas.surface_geometry[1])
    if atlas.surface_space is None or atlas.surface_density is None:
        return None
    return fetch_standard_surface_meshes(atlas.surface_space, atlas.surface_density, mesh_kind)


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


def surface_face_colors(face_values: np.ndarray, triangles: np.ndarray, coords: np.ndarray, *, cmap, norm) -> np.ndarray:
    """Return shaded RGBA face colors for one array of triangle values."""

    colors = np.asarray(cmap(norm(np.nan_to_num(face_values, nan=0.0))), dtype=float)
    shading = triangle_shading(coords, triangles)
    colors[:, :3] = np.clip(colors[:, :3] * shading[:, None], 0.0, 1.0)
    colors[~np.isfinite(face_values)] = np.array([0.92, 0.94, 0.97, 1.0], dtype=float)
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


def triangle_shading(coords: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Return subtle Lambertian-like shading values for each surface triangle."""

    tri_coords = coords[triangles]
    vec1 = tri_coords[:, 1] - tri_coords[:, 0]
    vec2 = tri_coords[:, 2] - tri_coords[:, 0]
    normals = np.cross(vec1, vec2)
    lengths = np.linalg.norm(normals, axis=1, keepdims=True)
    normals = np.divide(
        normals,
        lengths,
        out=np.zeros_like(normals),
        where=lengths > 0,
    )
    light = np.array([0.35, -0.25, 0.90], dtype=float)
    light /= np.linalg.norm(light)
    lambert = np.clip(normals @ light, -1.0, 1.0)
    return 0.93 + 0.08 * ((lambert + 1.0) / 2.0)


def surface_view(
    ax,
    coords: np.ndarray,
    triangles: np.ndarray,
    vertex_values: np.ndarray,
    *,
    azim: float,
    cmap,
    norm,
    elev: float = 8.0,
) -> None:
    """Render one cortical surface panel."""

    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    face_values = triangle_face_values(vertex_values, triangles)
    face_colors = surface_face_colors(face_values, triangles, coords, cmap=cmap, norm=norm)
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
    ax.view_init(elev=elev, azim=azim)
    ax.set_proj_type("ortho")
    ax.set_axis_off()
    extent = np.ptp(coords, axis=0)
    aspect = np.maximum(extent, 1e-3).astype(float)
    aspect[2] *= 1.08
    ax.set_box_aspect(tuple(aspect))
    center = coords.mean(axis=0)
    half_extent = (extent / 2.0) * np.array([1.06, 1.06, 1.12], dtype=float)
    ax.set_xlim(center[0] - half_extent[0], center[0] + half_extent[0])
    ax.set_ylim(center[1] - half_extent[1], center[1] + half_extent[1])
    ax.set_zlim(center[2] - half_extent[2], center[2] + half_extent[2])
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

    atlas = get_atlas(atlas_id)
    if atlas.surface_paths is None:
        return None

    hemi_frames = surface_value_frames(table)
    if not hemi_frames:
        return None

    inflated_paths = surface_mesh_paths(atlas_id)
    if inflated_paths is None:
        return None

    return _plot_cortical_surface_map_matplotlib(
        atlas=atlas,
        hemi_frames=hemi_frames,
        inflated_paths=inflated_paths,
        table=table,
        value_column=value_column,
        title=title,
        output_path=output_path,
    )


def plot_cortical_surface_map_brainspace(
    table: pd.DataFrame,
    *,
    atlas_id: str,
    value_column: str,
    title: str,
    output_path: Path,
) -> Path | None:
    """Render one regional table with BrainSpace when the backend is available."""

    if importlib.util.find_spec("brainspace") is None:
        return None

    atlas = get_atlas(atlas_id)
    if atlas.surface_paths is None:
        return None

    hemi_frames = surface_value_frames(table)
    if not hemi_frames:
        return None

    mesh_paths = surface_mesh_paths(atlas_id, mesh_kind="pial")
    if mesh_paths is None:
        return None

    try:
        vertex_payload: dict[str, np.ndarray] = {}
        if "left" in hemi_frames:
            label_lh, names_lh = load_surface_parcellation(str(atlas.surface_paths[0]))
            vertex_payload["left"] = vertex_values_for_hemisphere(
                hemi_frames["left"],
                value_column=value_column,
                label_array=label_lh,
                code_to_name=names_lh,
            )
        if "right" in hemi_frames and atlas.surface_paths[1] is not None:
            label_rh, names_rh = load_surface_parcellation(str(atlas.surface_paths[1]))
            vertex_payload["right"] = vertex_values_for_hemisphere(
                hemi_frames["right"],
                value_column=value_column,
                label_array=label_rh,
                code_to_name=names_rh,
            )
    except Exception as exc:  # pragma: no cover - fallback safety
        logger.warning("BrainSpace cortical rendering setup failed: %s", exc)
        return None

    if not vertex_payload:
        return None

    finite_values = table[value_column].to_numpy(dtype=float, copy=False)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size == 0:
        return None
    vmax = _display_vmax(finite_values, value_column=value_column)

    full_panel_labels: list[str] = []
    if "left" in vertex_payload:
        full_panel_labels.extend(["Left lateral", "Left medial"])
    if "right" in vertex_payload:
        full_panel_labels.extend(["Right medial", "Right lateral"])
    with tempfile.TemporaryDirectory(prefix="imt-brainspace-") as tmpdir:
        payload_path = Path(tmpdir) / "payload.npz"
        raw_output_path = Path(tmpdir) / "brainspace_raw.png"
        np.savez(
            payload_path,
            surf_lh_path=np.array(mesh_paths[0] if "left" in vertex_payload else "", dtype=object),
            surf_rh_path=np.array(mesh_paths[1] if "right" in vertex_payload else "", dtype=object),
            output_path=np.array(str(raw_output_path), dtype=object),
            vmax=np.array(float(vmax), dtype=float),
            has_left=np.array("left" in vertex_payload, dtype=bool),
            has_right=np.array("right" in vertex_payload, dtype=bool),
            left=vertex_payload.get("left", np.array([], dtype=float)),
            right=vertex_payload.get("right", np.array([], dtype=float)),
        )

        script = """
from pathlib import Path
import numpy as np
from brainspace.mesh.mesh_io import read_surface
from brainspace.plotting.surface_plotting import plot_surf
from brainspace.vtk_interface import serial_connect, wrap_vtk
from vtk import vtkPolyDataNormals

payload = np.load(r'__PAYLOAD__', allow_pickle=True)
output_path = Path(payload['output_path'].item())
vmax = float(payload['vmax'])
has_left = bool(payload['has_left'])
has_right = bool(payload['has_right'])
surf_lh_path = payload['surf_lh_path'].item()
surf_rh_path = payload['surf_rh_path'].item()

if not has_left and not has_right:
    raise RuntimeError('No cortical hemispheres available for BrainSpace rendering.')

surfs = {}
layout = []
views = []

def with_normals(surface):
    normals = wrap_vtk(vtkPolyDataNormals, splitting=False, featureAngle=0.1)
    return serial_connect(surface, normals)

if has_left:
    surf_lh = read_surface(surf_lh_path, return_data=True)
    surf_lh.append_array(payload['left'], name='imt_values', at='p')
    surf_lh = with_normals(surf_lh)
    surfs['lh'] = surf_lh
    layout.extend(['lh', 'lh'])
    views.extend(['lateral', 'medial'])

if has_right:
    surf_rh = read_surface(surf_rh_path, return_data=True)
    surf_rh.append_array(payload['right'], name='imt_values', at='p')
    surf_rh = with_normals(surf_rh)
    surfs['rh'] = surf_rh
    layout.extend(['rh', 'rh'])
    # In this custom row layout, BrainSpace's right-hemisphere camera names
    # render in the opposite visual order, so request lateral/medial here to
    # produce right-medial then right-lateral in the final figure.
    views.extend(['lateral', 'medial'])

plot_surf(
    surfs,
    [layout],
    array_name=[['imt_values'] * len(layout)],
    view=[views],
    color_bar=None,
    color_range=(-vmax, vmax),
    share='b',
    cmap='RdBu_r',
    nan_color=(0.92, 0.94, 0.97, 1.0),
    zoom=1.08,
    background=(1.0, 1.0, 1.0),
    size=(1450 if len(layout) == 4 else 960, 390),
    screenshot=True,
    filename=str(output_path),
    transparent_bg=False,
    scale=(2, 2),
    interactive=False,
    suppress_warnings=True,
    actor__ambient=0.16,
    actor__diffuse=0.86,
    actor__specular=0.03,
    actor__specularPower=6,
    actor__interpolation='phong',
)
"""
        script = script.replace("__PAYLOAD__", str(payload_path))
        env = os.environ.copy()
        env.setdefault("MPLCONFIGDIR", str(Path(tmpdir) / "mplconfig"))
        env.setdefault("PYTHONUNBUFFERED", "1")
        try:
            result = subprocess.run(
                [sys.executable, "-c", script],
                env=env,
                capture_output=True,
                text=True,
                timeout=90,
            )
        except subprocess.TimeoutExpired:
            logger.warning("BrainSpace cortical rendering timed out and was skipped.")
            return None

        if result.returncode != 0:
            stderr = result.stderr.strip().splitlines()
            detail = stderr[-1] if stderr else f"exit code {result.returncode}"
            logger.warning("BrainSpace cortical rendering failed and was skipped: %s", detail)
            return None

        if raw_output_path.exists():
            return _compose_brainspace_surface_figure(
                raw_output_path=raw_output_path,
                output_path=output_path,
                title=title,
                value_column=value_column,
                vmax=vmax,
                panel_labels=full_panel_labels,
            )
    return None


def _compose_brainspace_surface_figure(
    *,
    raw_output_path: Path,
    output_path: Path,
    title: str,
    value_column: str,
    vmax: float,
    panel_labels: list[str],
) -> Path:
    """Wrap a raw BrainSpace screenshot in the shared publication-style frame."""

    _, plt = matplotlib_backend()

    image = plt.imread(raw_output_path)
    fig_width = 12.4 if len(panel_labels) == 4 else 7.4
    fig = plt.figure(figsize=(fig_width, 4.6))
    fig.patch.set_facecolor("white")

    image_ax = fig.add_axes([0.02, 0.20, 0.96, 0.62])
    image_ax.imshow(image)
    image_ax.axis("off")

    label_centers = np.linspace(0.14, 0.86, len(panel_labels))
    for x, label in zip(label_centers, panel_labels, strict=False):
        fig.text(
            x,
            0.84,
            label,
            ha="center",
            va="bottom",
            fontsize=9.5,
            color="#334155",
        )

    fig.text(0.5, 0.94, title, ha="center", va="top", fontweight="bold", fontsize=12.5)

    norm = _surface_norm(plt, vmax)
    cmap = plt.get_cmap("RdBu_r")
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])
    colorbar_ax = fig.add_axes([0.36, 0.08, 0.28, 0.034])
    colorbar = fig.colorbar(mappable, cax=colorbar_ax, orientation="horizontal")
    colorbar.set_ticks([-vmax, 0.0, vmax])
    colorbar.ax.set_xticklabels(_surface_tick_labels(vmax))
    colorbar.ax.tick_params(labelsize=8.5, pad=1, length=2.5)
    colorbar.set_label(_surface_colorbar_label(value_column), fontsize=9.5, labelpad=3)
    colorbar.outline.set_linewidth(0.6)
    return save_figure(fig, output_path)


def _plot_cortical_surface_map_matplotlib(
    *,
    atlas,
    hemi_frames: dict[str, pd.DataFrame],
    inflated_paths: tuple[str, str],
    table: pd.DataFrame,
    value_column: str,
    title: str,
    output_path: Path,
) -> Path | None:
    """Render one cortical map with a publication-oriented Matplotlib layout."""

    _, plt = matplotlib_backend()

    panels: list[tuple[str, str, float]] = []
    if "left" in hemi_frames:
        panels.extend([("left", "Left lateral", 180), ("left", "Left medial", 0)])
    if "right" in hemi_frames:
        panels.extend([("right", "Right medial", 180), ("right", "Right lateral", 0)])
    if not panels:
        return None

    finite_values = table[value_column].to_numpy(dtype=float, copy=False)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size == 0:
        return None
    vmax = _display_vmax(finite_values, value_column=value_column)

    norm = _surface_norm(plt, vmax)
    cmap = plt.get_cmap("RdBu_r")
    figure_width = 12.4 if len(panels) == 4 else 7.4
    fig = plt.figure(figsize=(figure_width, 3.8))
    fig.patch.set_facecolor("white")

    hemi_index = {"left": 0, "right": 1}
    positions = _surface_panel_positions(len(panels))
    for (hemi, panel_title, azim), position in zip(panels, positions, strict=False):
        ax = fig.add_axes(position, projection="3d")
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
        surface_view(ax, coords, triangles, vertex_values, azim=azim, cmap=cmap, norm=norm)
        fig.text(
            position[0] + position[2] / 2.0,
            0.855,
            panel_title,
            ha="center",
            va="bottom",
            fontsize=9.0,
            color="#334155",
        )

    fig.text(0.035, 0.935, title, ha="left", va="top", fontweight="bold", fontsize=12.5)
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])
    colorbar_ax = fig.add_axes([0.385, 0.085, 0.23, 0.036])
    colorbar = fig.colorbar(mappable, cax=colorbar_ax, orientation="horizontal")
    colorbar.set_ticks([-vmax, 0.0, vmax])
    colorbar.ax.set_xticklabels(_surface_tick_labels(vmax))
    colorbar.ax.tick_params(labelsize=8.5, pad=1, length=2.5)
    colorbar.set_label(_surface_colorbar_label(value_column), fontsize=9.5, labelpad=3)
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.outline.set_linewidth(0.6)
    return save_figure(fig, output_path)


def _surface_panel_positions(n_panels: int) -> list[list[float]]:
    """Return tightly spaced figure coordinates for the cortical surface panels."""

    if n_panels == 4:
        return [
            [0.015, 0.14, 0.285, 0.69],
            [0.245, 0.14, 0.285, 0.69],
            [0.47, 0.14, 0.285, 0.69],
            [0.695, 0.14, 0.285, 0.69],
        ]
    if n_panels == 2:
        return [
            [0.06, 0.14, 0.40, 0.69],
            [0.46, 0.14, 0.40, 0.69],
        ]
    width = 0.84 / max(n_panels, 1)
    return [[0.08 + index * width, 0.18, width, 0.62] for index in range(n_panels)]


def _surface_colorbar_label(value_column: str) -> str:
    """Return a human-readable colorbar label for cortical surface plots."""

    if value_column.endswith("_z"):
        return "Z score"
    return str(value_column).replace("_", " ")


def _surface_tick_labels(vmax: float) -> list[str]:
    """Return compact symmetric tick labels for the cortical surface colorbar."""

    return [_format_surface_tick_value(-vmax), "0", _format_surface_tick_value(vmax)]


def _format_surface_tick_value(value: float) -> str:
    """Format one surface colorbar tick, trimming redundant decimals."""

    if np.isclose(value, round(value)):
        return str(int(round(value)))
    decimals = 1 if abs(value) >= 2 else 2
    return f"{value:.{decimals}f}"


def _display_vmax(finite_values: np.ndarray, *, value_column: str) -> float:
    """Return the symmetric display limit used by brain and cortical maps."""

    vmax = float(np.nanpercentile(np.abs(finite_values), 98))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = float(np.nanmax(np.abs(finite_values))) if finite_values.size else 1.0
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    if value_column.endswith("_z") and vmax >= 1.0:
        vmax = max(float(np.floor(vmax)), 1.0)
    return vmax


def _surface_norm(plt, vmax: float):
    """Create the diverging normalization used by cortical surface plots."""

    from matplotlib.colors import TwoSlopeNorm

    return TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)


def atlas_volume_data(atlas_id: str) -> np.ndarray | None:
    """Load one packaged atlas volume as integer parcel codes."""

    atlas_path = atlas_volume_path(atlas_id)
    if atlas_path is None:
        return None
    return np.asarray(nib.load(atlas_path).get_fdata(), dtype=np.int32)


@lru_cache(maxsize=None)
def atlas_label_table(atlas_id: str) -> pd.DataFrame | None:
    """Load the packaged atlas label table for slice-selection heuristics."""

    atlas = get_atlas(atlas_id)
    if atlas.labels_path is None:
        return None
    labels = pd.read_csv(atlas.labels_path)
    if "Unnamed: 0" in labels.columns:
        labels = labels.drop(columns=["Unnamed: 0"])
    return labels


def brain_volume_from_table(table: pd.DataFrame, atlas_id: str, value_column: str) -> np.ndarray | None:
    """Project one regional table into a volumetric atlas array."""

    if "id" not in table.columns or value_column not in table.columns:
        return None
    atlas_data = atlas_volume_data(atlas_id)
    if atlas_data is None:
        return None
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


def slice_index_from_mask(mask: np.ndarray, axis: int) -> int:
    """Choose one slice index from a boolean atlas-support mask."""

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


def _render_brain_slice(
    ax,
    atlas_slice: np.ndarray,
    data_slice: np.ndarray,
    *,
    cmap,
    norm,
) -> None:
    """Render one slice with smooth parcel-like fills inspired by ggseg plots."""

    from scipy.ndimage import gaussian_filter, zoom

    scale = 6
    sigma = 1.2
    outline_color = "#cbd5e1"
    boundary_color = "#ffffff"

    support = zoom((atlas_slice > 0).astype(float), scale, order=3)
    support = gaussian_filter(support, sigma=sigma)
    support_max = float(np.nanmax(support))
    if support_max > 0:
        ax.contourf(
            support,
            levels=[0.5, support_max + 1e-6],
            colors=[BACKGROUND],
            antialiased=True,
        )

    region_ids = [int(region_id) for region_id in np.unique(atlas_slice) if int(region_id) > 0]
    for region_id in region_ids:
        region_mask = atlas_slice == region_id
        region_values = data_slice[region_mask]
        region_value = (
            float(np.nanmean(region_values))
            if np.any(np.isfinite(region_values))
            else np.nan
        )
        fill_color = MISSING if not np.isfinite(region_value) else cmap(norm(float(region_value)))

        region = zoom(region_mask.astype(float), scale, order=3)
        region = gaussian_filter(region, sigma=sigma)
        region_max = float(np.nanmax(region))
        if region_max <= 0:
            continue
        ax.contourf(
            region,
            levels=[0.5, region_max + 1e-6],
            colors=[fill_color],
            antialiased=True,
        )
        ax.contour(
            region,
            levels=[0.5],
            colors=[boundary_color],
            linewidths=0.35,
            alpha=0.95,
        )

    if support_max > 0:
        ax.contour(
            support,
            levels=[0.5],
            colors=[outline_color],
            linewidths=0.9,
        )
    ax.set_aspect("equal")
    ax.invert_yaxis()


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
    atlas_data = atlas_volume_data(atlas_id)
    if atlas_data is None:
        return None

    finite_values = volume[np.isfinite(volume)]
    if finite_values.size == 0:
        return None
    vmax = _display_vmax(finite_values, value_column=value_column)

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import TwoSlopeNorm

    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    cmap = plt.get_cmap("RdBu_r")
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.2))
    fig.patch.set_facecolor("white")
    atlas_mask_volume = atlas_data > 0
    axial_mask = atlas_mask_volume
    labels = atlas_label_table(atlas_id)
    if labels is not None and "structure" in labels.columns and "id" in labels.columns:
        subcortical_ids = labels.loc[
            labels["structure"].astype(str) != "cortex",
            "id",
        ].to_numpy(dtype=np.int32, copy=False)
        if subcortical_ids.size:
            subcortical_mask = np.isin(atlas_data, subcortical_ids)
            if np.any(subcortical_mask):
                axial_mask = subcortical_mask
    slice_specs = [
        ("Sagittal", 0, slice_index_from_mask(atlas_mask_volume, 0)),
        ("Coronal", 1, slice_index_from_mask(atlas_mask_volume, 1)),
        ("Axial", 2, slice_index_from_mask(axial_mask, 2)),
    ]
    for ax, (label, axis, index) in zip(axes, slice_specs, strict=False):
        atlas_slice = slice_view(atlas_data.astype(float), axis, index)
        data_slice = slice_view(volume, axis, index)
        _render_brain_slice(ax, atlas_slice, data_slice, cmap=cmap, norm=norm)
        ax.set_title(label, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_facecolor(BACKGROUND)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle(title, x=0.06, y=0.98, ha="left", fontweight="bold", fontsize=12)
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axes, pad=0.02, shrink=0.85)
    colorbar.set_label(_surface_colorbar_label(value_column))
    return save_figure(fig, output_path)
