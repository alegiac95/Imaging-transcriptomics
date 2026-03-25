from __future__ import annotations

import warnings
from functools import lru_cache

import numpy as np

from .config import DEFAULT_NULL_METHOD, DEFAULT_SEED, VALID_NULL_METHODS
from .surfaces import load_surface_parcellation


NULL_METHODS = VALID_NULL_METHODS
SURFACE_NULL_METHODS = {"vasa", "alexander_bloch", "moran"}


def _import_neuromaps_nulls():
    try:
        from neuromaps import nulls
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "neuromaps is required for surface-based spatial null models. "
            "Install imaging-transcriptomics[maps]."
        ) from exc
    return nulls


def _standardize_vector(values: np.ndarray) -> np.ndarray:
    vector = np.asarray(values, dtype=float).reshape(-1)
    centered = vector - np.nanmean(vector)
    scale = np.nanstd(centered, ddof=1)
    if not np.isfinite(scale) or scale == 0:
        return centered
    return centered / scale


def _group_index_arrays(groups: np.ndarray) -> tuple[np.ndarray, ...]:
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
    null_maps = np.asarray(null_maps, dtype=float)
    if null_maps.shape == (n_values, n_permutations):
        return null_maps
    if null_maps.shape == (n_permutations, n_values):
        return null_maps.T
    raise ValueError(
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
    if method not in SURFACE_NULL_METHODS:
        raise ValueError(f"Unsupported surface null method: {method}")
    atlas = selection.atlas
    if atlas.surface_paths is None:
        raise FileNotFoundError(
            f"Atlas '{atlas.id}' does not have surface parcellation files for cortical null models."
        )
    if atlas.surface_space is None or atlas.surface_density is None:
        raise ValueError(
            f"Atlas '{atlas.id}' is missing surface space/density metadata required for cortical null models."
        )
    if selection.hemisphere == "both" and atlas.surface_paths[1] is None:
        raise FileNotFoundError(
            f"Atlas '{atlas.id}' does not have a right-hemisphere surface parcellation file for bilateral null models."
        )

    nulls = _import_neuromaps_nulls()
    parcellation = _surface_parcellation(
        atlas.id,
        str(atlas.surface_paths[0]),
        str(atlas.surface_paths[1]),
        selection.hemisphere,
    )
    kwargs = dict(
        data=np.asarray(cortical_values, dtype=float),
        atlas=atlas.surface_space,
        density=atlas.surface_density,
        parcellation=parcellation,
        n_perm=n_permutations,
        seed=seed,
    )
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
    if null_method not in NULL_METHODS:
        valid = ", ".join(sorted(NULL_METHODS))
        raise ValueError(f"Unknown null method '{null_method}'. Expected one of: {valid}.")

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
                try:
                    permuted[cortical_idx, :] = generate_surface_nulls(
                        values[cortical_idx],
                        extracted.selection,
                        method=method,
                        n_permutations=n_permutations,
                        seed=seed,
                    )
                    resolved_method = method
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
                    raise RuntimeError(
                        f"Unable to generate cortical nulls with method '{null_method}'."
                    ) from last_error

    non_cortical_idx = np.flatnonzero(label_structures != "cortex")
    if non_cortical_idx.size:
        non_cortical_hemi = labels.iloc[non_cortical_idx]["hemisphere"].astype(str).to_numpy()
        permuted[non_cortical_idx, :] = shuffle_within_groups(
            values[non_cortical_idx],
            non_cortical_hemi,
            n_permutations,
            rng=rng,
        )
    return permuted, resolved_method
