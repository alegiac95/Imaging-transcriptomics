from __future__ import annotations

from pathlib import Path

from ._compat import suppress_pkg_resources_deprecation


def _import_neuromaps_images():
    try:
        with suppress_pkg_resources_deprecation():
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


def surface_geometry_paths(atlas, hemisphere: str) -> tuple[str, ...] | None:
    geometry = surface_geometry(atlas)
    if geometry is None:
        return None
    if hemisphere == "left":
        return (geometry[0],)
    return geometry


def load_surface_parcellation(atlas, hemisphere: str):
    paths = surface_paths(atlas, hemisphere)
    if all(Path(path).suffix == ".annot" for path in paths):
        annot_to_gifti = _import_neuromaps_images()
        return annot_to_gifti(paths)
    return paths


def _parcellation_vertex_count(item) -> int:
    import nibabel as nib

    if isinstance(item, nib.GiftiImage):
        data = item.agg_data()
    else:
        path = Path(item)
        if path.suffix == ".annot":
            from nibabel.freesurfer import read_annot

            labels, _, _ = read_annot(str(path))
            return int(labels.shape[0])
        data = nib.load(str(path)).agg_data()

    if isinstance(data, tuple):
        data = data[0]
    return int(data.shape[0])


def infer_surface_density(atlas, hemisphere: str) -> str | None:
    space = getattr(atlas, "surface_space", None)
    default_density = getattr(atlas, "surface_density", None)
    parcellation = load_surface_parcellation(atlas, hemisphere)
    if not parcellation:
        return default_density

    counts = {_parcellation_vertex_count(item) for item in parcellation}
    if len(counts) != 1:
        return default_density
    count = counts.pop()

    mappings = {
        "fsaverage": {
            642: "1k",
            2562: "3k",
            10242: "10k",
            40962: "41k",
            163842: "164k",
        },
        "fslr": {
            32492: "32k",
            163842: "164k",
        },
        "civet": {
            40962: "41k",
        },
    }
    if space is None:
        return default_density
    return mappings.get(str(space).lower(), {}).get(count, default_density)
