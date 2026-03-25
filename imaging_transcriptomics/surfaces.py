from __future__ import annotations

from pathlib import Path


def _import_neuromaps_images():
    try:
        from neuromaps.images import annot_to_gifti
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "neuromaps is required for surface atlas support. "
            "Install imaging-transcriptomics[maps]."
        ) from exc
    return annot_to_gifti


def surface_paths(atlas, hemisphere: str) -> tuple[str, ...]:
    paths = atlas.surface_paths
    if paths is None:
        raise FileNotFoundError(f"Atlas '{atlas.id}' does not define surface parcellation files.")
    if hemisphere == "left":
        return (str(paths[0]),)
    return tuple(str(path) for path in paths)


def surface_geometry(atlas) -> tuple[str, str] | None:
    geometry = atlas.surface_geometry
    if geometry is None:
        return None
    return tuple(str(path) for path in geometry)


def load_surface_parcellation(atlas, hemisphere: str):
    paths = surface_paths(atlas, hemisphere)
    if all(Path(path).suffix == ".annot" for path in paths):
        annot_to_gifti = _import_neuromaps_images()
        return annot_to_gifti(paths)
    return paths
