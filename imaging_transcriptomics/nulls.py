from __future__ import annotations

import warnings
from functools import lru_cache

import numpy as np

from ._logging import get_logger
from ._compat import suppress_pkg_resources_deprecation
from .config import DEFAULT_NULL_METHOD, DEFAULT_SEED, VALID_NULL_METHODS
from .exceptions import AtlasAssetError, NullModelError
from .surfaces import infer_surface_density, load_surface_parcellation, surface_geometry_paths


NULL_METHODS = VALID_NULL_METHODS
SURFACE_NULL_METHODS = {"vasa", "alexander_bloch", "moran"}
logger = get_logger(__name__)


def _import_neuromaps_nulls():
    """Import neuromaps null generators lazily and under warning suppression."""

    try:
        with suppress_pkg_resources_deprecation():
            from neuromaps import nulls
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "neuromaps is required for surface-based spatial null models and "
            "is part of the standard imaging-transcriptomics install. "
            "Reinstall the package if it is missing from the current "
            "environment."
        ) from exc
    return nulls


def _standardize_vector(values: np.ndarray) -> np.ndarray:
    """Center and scale a vector while ignoring NaNs."""

    vector = np.asarray(values, dtype=float).reshape(-1)
    centered = vector - np.nanmean(vector)
    scale = np.nanstd(centered, ddof=1)
    if not np.isfinite(scale) or scale == 0:
        return centered
    return centered / scale


def _group_index_arrays(groups: np.ndarray) -> tuple[np.ndarray, ...]:
    """Return index arrays for each unique grouping label."""

    groups = np.asarray(groups)
    _, inverse = np.unique(groups.astype(str), return_inverse=True)
    return tuple(np.flatnonzero(inverse == idx) for idx in range(inverse.max() + 1))


def shuffle_within_groups(
    values: np.ndarray,
    groups: np.ndarray,
    n_permutations: int,
    *,
    rng: np.random.Generator,
) -> np.ndarray:
    """Permute values independently within predefined region groups."""

    values = np.asarray(values, dtype=float).reshape(-1)
    permuted = np.repeat(values[:, None], n_permutations, axis=1)
    for idx in _group_index_arrays(groups):
        if idx.size < 2:
            continue
        order = np.argsort(rng.random((idx.size, n_permutations)), axis=0)
        permuted[idx, :] = values[idx][order]
    return permuted


@lru_cache(maxsize=None)
def _surface_parcellation(
    atlas_id: str,
    lh_surface_path: str,
    rh_surface_path: str | None,
    hemisphere: str,
):
    """Load and cache the surface parcellation used by neuromaps null models."""

    del atlas_id

    class _AtlasProxy:
        def __init__(self):
            if rh_surface_path is None:
                self.surface_paths = (lh_surface_path, lh_surface_path)
            else:
                self.surface_paths = (lh_surface_path, rh_surface_path)
            self.id = "surface"

    return load_surface_parcellation(_AtlasProxy(), hemisphere)


def _normalize_null_maps(null_maps: np.ndarray, n_values: int, n_permutations: int, method: str) -> np.ndarray:
    """Normalize neuromaps null output orientation to values-by-permutation."""

    null_maps = np.asarray(null_maps, dtype=float)
    if null_maps.shape == (n_values, n_permutations):
        return null_maps
    if null_maps.shape == (n_permutations, n_values):
        return null_maps.T
    raise NullModelError(
        f"Null model '{method}' returned an unexpected shape {null_maps.shape}; "
        f"expected ({n_values}, {n_permutations})."
    )


def generate_surface_nulls(
    cortical_values: np.ndarray,
    selection,
    *,
    method: str,
    n_permutations: int,
    seed: int,
) -> np.ndarray:
    """Generate cortical spatial nulls with a requested neuromaps method."""

    if method not in SURFACE_NULL_METHODS:
        raise NullModelError(f"Unsupported surface null method: {method}")
    atlas = selection.atlas
    if atlas.surface_paths is None:
        raise AtlasAssetError(
            f"Atlas '{atlas.id}' does not have surface parcellation files for cortical null models."
        )
    if atlas.surface_space is None or atlas.surface_density is None:
        raise NullModelError(
            f"Atlas '{atlas.id}' is missing surface space/density metadata required for cortical null models."
        )
    if selection.hemisphere == "both" and atlas.surface_paths[1] is None:
        raise AtlasAssetError(
            f"Atlas '{atlas.id}' does not have a right-hemisphere surface parcellation file for bilateral null models."
        )

    logger.info(
        "Generating %d cortical spatial nulls with '%s' for atlas '%s'.",
        n_permutations,
        method,
        atlas.id,
    )
    nulls = _import_neuromaps_nulls()
    parcellation = _surface_parcellation(
        atlas.id,
        str(atlas.surface_paths[0]),
        str(atlas.surface_paths[1]),
        selection.hemisphere,
    )
    geometry = surface_geometry_paths(atlas, selection.hemisphere)
    kwargs = dict(
        data=np.asarray(cortical_values, dtype=float),
        atlas=atlas.surface_space,
        density=infer_surface_density(atlas, selection.hemisphere),
        parcellation=parcellation,
        n_perm=n_permutations,
        seed=seed,
    )
    if geometry is not None:
        kwargs["surfaces"] = geometry
    if method == "vasa":
        generated = nulls.vasa(**kwargs)
    elif method == "alexander_bloch":
        generated = nulls.alexander_bloch(**kwargs)
    else:
        generated = nulls.moran(**kwargs, n_proc=1)
    return _normalize_null_maps(generated, cortical_values.shape[0], n_permutations, method)


def permute_scan_values(
    extracted,
    n_permutations: int,
    *,
    null_method: str = DEFAULT_NULL_METHOD,
    seed: int = DEFAULT_SEED,
) -> tuple[np.ndarray, str]:
    """Generate regional null maps for one extracted imaging vector.

    Cortical regions use the requested surface null model when possible, while
    non-cortical regions are permuted within hemisphere groups.
    """

    if null_method not in NULL_METHODS:
        valid = ", ".join(sorted(NULL_METHODS))
        raise NullModelError(f"Unknown null method '{null_method}'. Expected one of: {valid}.")

    logger.info(
        "Generating %d permuted maps with null method '%s'.",
        n_permutations,
        null_method,
    )
    values = _standardize_vector(extracted.values)
    labels = extracted.labels.reset_index(drop=True)
    permuted = np.zeros((values.shape[0], n_permutations), dtype=float)
    rng = np.random.default_rng(seed)
    resolved_method = null_method

    label_structures = labels["structure"].to_numpy(dtype=object, copy=False)
    cortical_idx = np.flatnonzero(label_structures == "cortex")
    if cortical_idx.size:
        cortical_hemi = labels.iloc[cortical_idx]["hemisphere"].astype(str).to_numpy()
        if null_method == "random":
            logger.info(
                "Using within-hemisphere random shuffling for %d cortical regions.",
                cortical_idx.size,
            )
            permuted[cortical_idx, :] = shuffle_within_groups(
                values[cortical_idx],
                cortical_hemi,
                n_permutations,
                rng=rng,
            )
        else:
            method_order = [null_method] if null_method != "auto" else ["vasa", "alexander_bloch"]
            last_error = None
            for method in method_order:
                logger.info(
                    "Trying cortical null method '%s' for %d cortical regions.",
                    method,
                    cortical_idx.size,
                )
                try:
                    permuted[cortical_idx, :] = generate_surface_nulls(
                        values[cortical_idx],
                        extracted.selection,
                        method=method,
                        n_permutations=n_permutations,
                        seed=seed,
                    )
                    resolved_method = method
                    logger.info("Using cortical null method '%s'.", method)
                    break
                except Exception as exc:
                    last_error = exc
            else:
                if null_method == "auto":
                    warnings.warn(
                        "Falling back to within-hemisphere random cortical permutations because surface nulls are unavailable: "
                        f"{last_error}",
                        RuntimeWarning,
                    )
                    permuted[cortical_idx, :] = shuffle_within_groups(
                        values[cortical_idx],
                        cortical_hemi,
                        n_permutations,
                        rng=rng,
                    )
                    resolved_method = "random"
                else:
                    raise NullModelError(
                        f"Unable to generate cortical nulls with method '{null_method}'."
                    ) from last_error

    non_cortical_idx = np.flatnonzero(label_structures != "cortex")
    if non_cortical_idx.size:
        non_cortical_hemi = labels.iloc[non_cortical_idx]["hemisphere"].astype(str).to_numpy()
        logger.info(
            "Generating within-hemisphere permutations for %d non-cortical regions.",
            non_cortical_idx.size,
        )
        permuted[non_cortical_idx, :] = shuffle_within_groups(
            values[non_cortical_idx],
            non_cortical_hemi,
            n_permutations,
            rng=rng,
        )
    logger.info("Finished generating permuted maps with resolved method '%s'.", resolved_method)
    return permuted, resolved_method
