from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .atlas_registry import normalize_atlas_id
from .exceptions import ConfigurationError
from .models import AnalysisMethod, HemisphereMode, NullMethod, RegionScope
from .validation import ensure_choice, ensure_positive_int, ensure_probability


DEFAULT_PERMUTATIONS = 1000
DEFAULT_NULL_METHOD = "auto"
DEFAULT_SEED = 1234
DEFAULT_N_JOBS = 1
DEFAULT_GENESET_ORGANISM = "Human"
DEFAULT_ENRICHMENT_METHOD = "ensemble"
VALID_HEMISPHERES = {"left", "both"}
VALID_REGIONS = {"default", "all", "cort", "cort+sub"}
VALID_METHODS = {"corr", "pls"}
VALID_NULL_METHODS = {"auto", "vasa", "alexander_bloch", "moran", "random"}
VALID_ENRICHMENT_METHODS = {"ensemble", "gsea", "ora", "none"}


@dataclass(frozen=True)
class RunConfig:
    """Validated v2 run configuration shared across the API and CLI."""

    method: AnalysisMethod
    atlas: str = "dk"
    hemisphere: HemisphereMode = "left"
    regions: RegionScope = "default"
    source_space: str | None = None
    n_permutations: int = DEFAULT_PERMUTATIONS
    null_method: NullMethod = DEFAULT_NULL_METHOD
    output_dir: Path | None = None
    enrichment_method: str = DEFAULT_ENRICHMENT_METHOD
    run_gsea: bool = False
    gene_set: str = "lake"
    geneset_organism: str = DEFAULT_GENESET_ORGANISM
    ora_p_threshold: float | None = None
    n_components: int | None = None
    var: float | None = None
    seed: int = DEFAULT_SEED
    n_jobs: int = DEFAULT_N_JOBS


def ensure_output_dir(output_dir: str | Path | None) -> Path | None:
    """Create an output directory if requested and return it as a Path."""

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
    regions: str = "default",
    source_space: str | None = None,
    n_permutations: int = DEFAULT_PERMUTATIONS,
    null_method: str = DEFAULT_NULL_METHOD,
    output_dir: str | Path | None = None,
    enrichment_method: str | None = None,
    run_gsea: bool = False,
    gene_set: str = "lake",
    geneset_organism: str = DEFAULT_GENESET_ORGANISM,
    ora_p_threshold: float | None = None,
    n_components: int | None = None,
    var: float | None = None,
    seed: int = DEFAULT_SEED,
    n_jobs: int = DEFAULT_N_JOBS,
) -> RunConfig:
    """Validate common analysis options and return a normalized run config.

    This helper is shared by the Python API and the CLI so both entry points use
    the same validation rules for atlas choice, region scope, permutation count,
    null model, enrichment settings, and PLS-specific options.
    """

    method = ensure_choice(method, VALID_METHODS, name="method")
    atlas_id = normalize_atlas_id(atlas)
    hemisphere = ensure_choice(hemisphere, VALID_HEMISPHERES, name="hemisphere")
    regions = ensure_choice(regions, VALID_REGIONS, name="regions")
    n_permutations = ensure_positive_int(n_permutations, name="n_permutations")
    null_method = ensure_choice(null_method, VALID_NULL_METHODS, name="null_method")
    ora_p_threshold = None if ora_p_threshold is None else ensure_probability(ora_p_threshold, name="ora_p_threshold")

    seed = int(seed)
    n_jobs = ensure_positive_int(n_jobs, name="n_jobs")
    resolved_output = ensure_output_dir(output_dir)

    if enrichment_method is not None:
        enrichment_method = ensure_choice(enrichment_method, VALID_ENRICHMENT_METHODS, name="enrichment_method")
    elif run_gsea:
        enrichment_method = "gsea"
    elif ora_p_threshold is not None:
        enrichment_method = "ora"
    else:
        enrichment_method = DEFAULT_ENRICHMENT_METHOD

    if enrichment_method == "ora" and ora_p_threshold is None:
        ora_p_threshold = 0.05

    if method == "pls":
        if n_components is None and var is None:
            raise ConfigurationError("PLS runs require either n_components or var.")
        if n_components is not None:
            n_components = ensure_positive_int(n_components, name="n_components")
        if var is not None:
            var = ensure_probability(var, name="var")

    return RunConfig(
        method=method,
        atlas=atlas_id,
        hemisphere=hemisphere,
        regions=regions,
        source_space=source_space,
        n_permutations=n_permutations,
        null_method=null_method,
        output_dir=resolved_output,
        enrichment_method=enrichment_method,
        run_gsea=bool(run_gsea),
        gene_set=gene_set,
        geneset_organism=str(geneset_organism),
        ora_p_threshold=ora_p_threshold,
        n_components=n_components,
        var=var,
        seed=seed,
        n_jobs=n_jobs,
    )
