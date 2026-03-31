from __future__ import annotations

import argparse
import os
import sys
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


sys.path.insert(0, str(_repo_root()))

from imaging_transcriptomics.atlas_registry import get_atlas, list_atlases  # noqa: E402
from imaging_transcriptomics.outputs.brain import (  # noqa: E402
    load_surface_mesh,
    load_surface_parcellation,
    surface_mesh_paths,
    surface_value_frames,
    surface_view,
    vertex_values_for_hemisphere,
)
from imaging_transcriptomics.outputs.common import matplotlib_backend  # noqa: E402


def gallery_mesh_paths(atlas_id: str, *, mesh_kind: str = "inflated") -> tuple[str, str] | None:
    """Return a consistent reference surface mesh for atlas gallery renders."""

    atlas = get_atlas(atlas_id)
    if atlas.surface_space is not None and atlas.surface_density is not None:
        cache_dir = _repo_root() / ".cache" / "neuromaps"
        cache_dir.mkdir(parents=True, exist_ok=True)
        try:
            from neuromaps.datasets import fetch_atlas
            atlas_files = fetch_atlas(
                atlas.surface_space,
                atlas.surface_density,
                data_dir=str(cache_dir),
                verbose=0,
            )
            candidates: list[object] = []
            if mesh_kind == "inflated":
                candidates.extend(
                    [
                        getattr(atlas_files, "inflated", None),
                        getattr(atlas_files, "pial", None),
                        getattr(atlas_files, "midthickness", None),
                        getattr(atlas_files, "white", None),
                    ]
                )
            else:
                candidates.extend(
                    [
                        getattr(atlas_files, "pial", None),
                        getattr(atlas_files, "inflated", None),
                        getattr(atlas_files, "midthickness", None),
                        getattr(atlas_files, "white", None),
                    ]
                )
            mesh = next((candidate for candidate in candidates if candidate is not None), None)
            if mesh is not None:
                return str(mesh[0]), str(mesh[1])
        except Exception:
            pass
    return surface_mesh_paths(atlas_id, mesh_kind=mesh_kind)


def destrieux_surface_parcellation() -> tuple[np.ndarray, dict[int, str]]:
    """Load the standard fsaverage5 Destrieux surface parcellation."""

    from nilearn.datasets import fetch_atlas_surf_destrieux

    cache_dir = _repo_root() / ".cache" / "nilearn"
    cache_dir.mkdir(parents=True, exist_ok=True)
    atlas = fetch_atlas_surf_destrieux(data_dir=str(cache_dir), verbose=0)
    labels = np.asarray(atlas.map_left, dtype=np.int32)
    code_to_name = {int(code): str(name) for code, name in enumerate(atlas.labels)}
    return labels, code_to_name


def load_labels(atlas_id: str) -> pd.DataFrame:
    """Load one packaged atlas label table with a consistent schema."""

    atlas = get_atlas(atlas_id)
    if atlas.labels_path is None:
        raise ValueError(f"Atlas '{atlas_id}' does not provide a packaged label table.")
    labels = pd.read_csv(atlas.labels_path)
    if "Unnamed: 0" in labels.columns:
        labels = labels.drop(columns=["Unnamed: 0"])
    return labels


def cortical_gallery_table(atlas_id: str) -> pd.DataFrame:
    """Return cortical atlas rows with a synthetic display value per parcel."""

    labels = load_labels(atlas_id)
    cortical = labels.loc[labels["structure"].astype(str) == "cortex"].copy()
    if cortical.empty:
        raise ValueError(f"Atlas '{atlas_id}' does not expose cortical parcels for gallery rendering.")

    # Spread parcel ids smoothly across the hue range inside each hemisphere so
    # both sides display a comparable amount of color variation without
    # implying any biological ordering.
    cortical["gallery_value"] = np.nan
    for hemisphere in ("L", "R"):
        mask = cortical["hemisphere"].astype(str) == hemisphere
        n_rows = int(mask.sum())
        if n_rows == 0:
            continue
        cortical.loc[mask, "gallery_value"] = np.linspace(0.0, 1.0, n_rows, endpoint=True)
    return cortical


def _crop_white_margin(image: np.ndarray, *, pad: int = 10) -> np.ndarray:
    """Crop a mostly white screenshot down to its visible surface content."""

    if image.ndim == 2:
        mask = image < 0.995
    else:
        rgb = image[..., :3]
        alpha = image[..., 3] if image.shape[-1] == 4 else np.ones(image.shape[:2], dtype=float)
        mask = (alpha > 0.01) & np.any(rgb < 0.995, axis=-1)
    if not np.any(mask):
        return image

    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    row_start = max(int(rows[0]) - pad, 0)
    row_stop = min(int(rows[-1]) + pad + 1, image.shape[0])
    col_start = max(int(cols[0]) - pad, 0)
    col_stop = min(int(cols[-1]) + pad + 1, image.shape[1])
    return image[row_start:row_stop, col_start:col_stop]


def _render_brainspace_left_lateral(
    *,
    atlas_id: str,
    vertex_values: np.ndarray,
    output_path: Path,
) -> Path | None:
    """Render a single left-lateral atlas preview with BrainSpace."""

    mesh_paths = gallery_mesh_paths(atlas_id, mesh_kind="inflated")
    if mesh_paths is None:
        return None

    with tempfile.TemporaryDirectory(prefix="imt-atlas-gallery-") as tmpdir:
        payload_path = Path(tmpdir) / "payload.npz"
        raw_output_path = Path(tmpdir) / "brainspace_left.png"
        np.savez(
            payload_path,
            surf_lh_path=np.array(mesh_paths[0], dtype=object),
            output_path=np.array(str(raw_output_path), dtype=object),
            left=vertex_values,
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
surf_lh_path = payload['surf_lh_path'].item()

surf_lh = read_surface(surf_lh_path, return_data=True)
surf_lh.append_array(payload['left'], name='atlas_preview', at='p')
normals = wrap_vtk(vtkPolyDataNormals, splitting=False, featureAngle=0.1)
surf_lh = serial_connect(surf_lh, normals)

    plot_surf(
        {'lh': surf_lh},
        [['lh']],
    array_name=[['atlas_preview']],
    view=[['lateral']],
    color_bar=None,
    color_range=(0.0, 1.0),
    share='b',
    cmap='viridis',
    nan_color=(1.0, 1.0, 1.0, 1.0),
    zoom=1.14,
    background=(1.0, 1.0, 1.0),
    size=(520, 380),
    screenshot=True,
    filename=str(output_path),
    transparent_bg=True,
    scale=(2, 2),
    interactive=False,
    suppress_warnings=True,
    actor__ambient=0.28,
    actor__diffuse=0.72,
    actor__specular=0.01,
    actor__specularPower=3,
    actor__interpolation='phong',
)
"""
        script = script.replace("__PAYLOAD__", str(payload_path))
        env = os.environ.copy()
        env.setdefault("MPLCONFIGDIR", str(Path(tmpdir) / "mplconfig"))
        env.setdefault("PYTHONUNBUFFERED", "1")
        result = subprocess.run(
            [sys.executable, "-c", script],
            env=env,
            capture_output=True,
            text=True,
            timeout=90,
        )
        if result.returncode != 0 or not raw_output_path.exists():
            return None

        _, plt = matplotlib_backend()
        cropped = _crop_white_margin(plt.imread(raw_output_path), pad=12)
        plt.imsave(output_path, cropped)
        return output_path


def _render_matplotlib_left_lateral(
    *,
    atlas_id: str,
    vertex_values: np.ndarray,
    output_path: Path,
) -> Path | None:
    """Render a single left-lateral atlas preview with Matplotlib."""

    _, plt = matplotlib_backend()
    mesh_paths = gallery_mesh_paths(atlas_id, mesh_kind="inflated")
    if mesh_paths is None:
        return None

    from matplotlib.colors import Normalize

    coords, triangles = load_surface_mesh(mesh_paths[0])
    fig = plt.figure(figsize=(3.6, 2.8), facecolor=(1.0, 1.0, 1.0, 0.0))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection="3d")
    surface_view(
        ax,
        coords,
        triangles,
        vertex_values,
        azim=180.0,
        cmap=plt.get_cmap("viridis"),
        norm=Normalize(vmin=0.0, vmax=1.0),
        elev=10.0,
    )
    ax.set_facecolor((1.0, 1.0, 1.0, 0.0))
    tmp_output = output_path.with_suffix(".tmp.png")
    fig.savefig(tmp_output, dpi=220, transparent=True)
    plt.close(fig)
    cropped = _crop_white_margin(plt.imread(tmp_output), pad=10)
    plt.imsave(output_path, cropped)
    tmp_output.unlink(missing_ok=True)
    return output_path


def render_atlas_gallery_figure(atlas_id: str, output_path: Path) -> Path | None:
    """Render one atlas preview as a left-lateral cortical gallery figure."""

    atlas = get_atlas(atlas_id)
    if atlas.surface_paths is None and atlas_id != "destrieux":
        return None

    table = cortical_gallery_table(atlas_id)
    hemi_frames = surface_value_frames(table)
    if "left" not in hemi_frames:
        return None

    if atlas_id == "destrieux":
        label_array, code_to_name = destrieux_surface_parcellation()
    else:
        label_array, code_to_name = load_surface_parcellation(str(atlas.surface_paths[0]))
    vertex_values = vertex_values_for_hemisphere(
        hemi_frames["left"],
        value_column="gallery_value",
        label_array=label_array,
        code_to_name=code_to_name,
    )
    rendered = _render_brainspace_left_lateral(
        atlas_id=atlas_id,
        vertex_values=vertex_values,
        output_path=output_path,
    )
    if rendered is not None:
        return rendered
    return _render_matplotlib_left_lateral(
        atlas_id=atlas_id,
        vertex_values=vertex_values,
        output_path=output_path,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render atlas gallery figures for the docs site.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_repo_root() / "docs" / "chapters" / "images" / "atlases",
        help="Directory where the rendered atlas gallery PNGs should be written.",
    )
    parser.add_argument(
        "--atlases",
        nargs="*",
        default=[atlas.id for atlas in list_atlases(packaged_only=True)],
        help="Atlas ids to render. Defaults to all packaged atlases.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for atlas_id in args.atlases:
        output_path = args.output_dir / f"{atlas_id}_gallery.png"
        rendered = render_atlas_gallery_figure(atlas_id, output_path)
        if rendered is None:
            print(f"skipped {atlas_id}: no cortical surface preview available")
        else:
            print(rendered)


if __name__ == "__main__":
    main()
