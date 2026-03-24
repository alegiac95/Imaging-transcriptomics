from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .atlas_registry import normalize_atlas_id
from .models import AnalysisMethod, HemisphereMode, RegionScope


DEFAULT_PERMUTATIONS = 1000
DEFAULT_NULL_METHOD = "auto"
DEFAULT_SEED = 1234
VALID_HEMISPHERES = {"left", "both"}
VALID_REGIONS = {"all", "cort", "cort+sub"}
VALID_METHODS = {"corr", "pls"}
VALID_NULL_METHODS = {"auto", "vasa", "alexander_bloch", "moran", "random"}


@dataclass(frozen=True)
class RunConfig:
    """Validated v2 run configuration shared across the API and CLI."""

    method: AnalysisMethod
    atlas: str = "dk"
    hemisphere: HemisphereMode = "left"
    regions: RegionScope = "all"
    source_space: str | None = None
    n_permutations: int = DEFAULT_PERMUTATIONS
    null_method: str = DEFAULT_NULL_METHOD
    output_dir: Path | None = None
    run_gsea: bool = False
    gene_set: str = "lake"
    n_components: int | None = None
    var: float | None = None
    seed: int = DEFAULT_SEED


def ensure_output_dir(output_dir: str | Path | None) -> Path | None:
    if output_dir is None:
        return None
    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)
    return path


def build_run_config(
    method: str,
    *,
    atlas: str = "dk",
    hemisphere: str = "left",
    regions: str = "all",
    source_space: str | None = None,
    n_permutations: int = DEFAULT_PERMUTATIONS,
    null_method: str = DEFAULT_NULL_METHOD,
    output_dir: str | Path | None = None,
    run_gsea: bool = False,
    gene_set: str = "lake",
    n_components: int | None = None,
    var: float | None = None,
    seed: int = DEFAULT_SEED,
) -> RunConfig:
    method = str(method).lower()
    if method not in VALID_METHODS:
        valid = ", ".join(sorted(VALID_METHODS))
        raise ValueError(f"Unknown method '{method}'. Expected one of: {valid}.")

    atlas_id = normalize_atlas_id(atlas)
    if hemisphere not in VALID_HEMISPHERES:
        raise ValueError("hemisphere must be either 'left' or 'both'.")
    if regions not in VALID_REGIONS:
        raise ValueError("regions must be one of 'all', 'cort', or 'cort+sub'.")

    n_permutations = int(n_permutations)
    if n_permutations < 1:
        raise ValueError("n_permutations must be at least 1.")
    if null_method not in VALID_NULL_METHODS:
        valid = ", ".join(sorted(VALID_NULL_METHODS))
        raise ValueError(f"Unknown null_method '{null_method}'. Expected one of: {valid}.")

    seed = int(seed)
    resolved_output = ensure_output_dir(output_dir)

    if method == "pls":
        if n_components is None and var is None:
            raise ValueError("PLS runs require either n_components or var.")
        if n_components is not None and int(n_components) < 1:
            raise ValueError("n_components must be at least 1.")
        if var is not None and not 0 < float(var) <= 1:
            raise ValueError("var must be in the interval (0, 1].")

    return RunConfig(
        method=method,
        atlas=atlas_id,
        hemisphere=hemisphere,
        regions=regions,
        source_space=source_space,
        n_permutations=n_permutations,
        null_method=null_method,
        output_dir=resolved_output,
        run_gsea=bool(run_gsea),
        gene_set=gene_set,
        n_components=None if n_components is None else int(n_components),
        var=None if var is None else float(var),
        seed=seed,
    )
