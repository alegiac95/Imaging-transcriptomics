"""Public functional API for imaging transcriptomics analyses."""

from __future__ import annotations

from .config import DEFAULT_NULL_METHOD, DEFAULT_PERMUTATIONS, RunConfig, build_run_config
from .gene_query import run_gene as _run_gene
from .models import CorrelationResult, PLSResult
from .nulls import permute_scan_values as _permute_scan_values
from .workflows.correlation import run_corr_configured as _workflow_run_corr_configured
from .workflows.pls import run_pls_configured as _workflow_run_pls_configured


def _run_corr_configured(data, config: RunConfig, *, input_rh=None) -> CorrelationResult:
    """Execute a configured correlation workflow using the public API seam."""

    return _workflow_run_corr_configured(
        data,
        config,
        input_rh=input_rh,
        permute_scan_values_fn=_permute_scan_values,
    )


def _run_pls_configured(data, config: RunConfig, *, input_rh=None) -> PLSResult:
    """Execute a configured PLS workflow using the public API seam."""

    return _workflow_run_pls_configured(
        data,
        config,
        input_rh=input_rh,
        permute_scan_values_fn=_permute_scan_values,
    )


def run_analysis(data, config: RunConfig, *, input_rh=None) -> CorrelationResult | PLSResult:
    """Run either correlation or PLS analysis from a validated configuration."""

    if config.method == "corr":
        return _run_corr_configured(data, config, input_rh=input_rh)
    return _run_pls_configured(data, config, input_rh=input_rh)


def run_corr(
    data,
    *,
    atlas: str = "dk",
    hemisphere: str = "left",
    regions: str = "default",
    source_space: str | None = None,
    input_rh=None,
    n_permutations: int = DEFAULT_PERMUTATIONS,
    null_method: str = DEFAULT_NULL_METHOD,
    output_dir=None,
    enrichment_method: str | None = None,
    run_gsea: bool = False,
    gene_set: str = "lake",
    geneset_organism: str = "Human",
    ora_p_threshold: float | None = None,
    seed: int = 1234,
    n_jobs: int = 1,
) -> CorrelationResult:
    """Run spatial correlation between an imaging map and atlas gene expression.

    Parameters
    ----------
    data
        Imaging input. This can be a regional vector, a text file containing one
        regional value per row, a volumetric NIfTI map, or a left/right surface
        pair that can be resampled with ``neuromaps``.
    atlas
        Atlas identifier registered in the packaged atlas registry.
    hemisphere
        Hemisphere subset to analyse. Use ``"left"`` for the historical
        left-only workflow or ``"both"`` for mirrored bilateral expression.
    regions
        Atlas region subset. ``"default"`` and ``"cort"`` keep cortex only;
        ``"all"`` and ``"cort+sub"`` include the packaged aseg add-on.
    source_space
        Declared space of the input image when it is not already on the packaged
        MNI atlas grid.
    input_rh
        Optional right-hemisphere surface file when ``data`` points to the left
        hemisphere file.
    n_permutations
        Number of spatial null permutations used for gene-wise statistics.
    null_method
        Spatial null model to use. ``"auto"`` chooses a supported method based
        on the atlas and available dependencies.
    output_dir
        Optional output directory for TSV, metadata, and plot files.
    enrichment_method
        Enrichment backend to apply after the gene-wise statistics are
        computed. Use ``"ensemble"`` for phenotype-null category enrichment,
        ``"gsea"`` for preranked GSEA, ``"ora"`` for over-representation
        analysis, or ``"none"`` to skip enrichment. When omitted, the default
        is ``"ensemble"`` unless the legacy ``run_gsea`` or
        ``ora_p_threshold`` arguments imply another choice.
    run_gsea
        Legacy compatibility flag that forces the preranked GSEA backend.
    gene_set
        Packaged geneset name, Enrichr library name, or GMT path used by GSEA
        and ORA.
    geneset_organism
        Organism used when ``gene_set`` refers to an Enrichr/GSEApy library
        name rather than a packaged or local GMT file.
    ora_p_threshold
        If provided, run ORA on genes whose raw permutation p-value is at or
        below this threshold, split into positive and negative gene lists.
    seed
        Random seed for the spatial null generation.
    n_jobs
        Reserved for API symmetry with PLS. Correlation currently uses a
        vectorized permutation path.

    Returns
    -------
    CorrelationResult
        Structured result object containing metadata, regional values, gene
        statistics, and optional enrichment outputs.
    """

    config = build_run_config(
        "corr",
        atlas=atlas,
        hemisphere=hemisphere,
        regions=regions,
        source_space=source_space,
        n_permutations=n_permutations,
        null_method=null_method,
        output_dir=output_dir,
        enrichment_method=enrichment_method,
        run_gsea=run_gsea,
        gene_set=gene_set,
        geneset_organism=geneset_organism,
        ora_p_threshold=ora_p_threshold,
        seed=seed,
        n_jobs=n_jobs,
    )
    return _run_corr_configured(data, config, input_rh=input_rh)


def run_pls(
    data,
    *,
    atlas: str = "dk",
    hemisphere: str = "left",
    regions: str = "default",
    source_space: str | None = None,
    input_rh=None,
    n_components: int | None = None,
    var: float | None = None,
    n_permutations: int = DEFAULT_PERMUTATIONS,
    null_method: str = DEFAULT_NULL_METHOD,
    output_dir=None,
    enrichment_method: str | None = None,
    run_gsea: bool = False,
    gene_set: str = "lake",
    geneset_organism: str = "Human",
    ora_p_threshold: float | None = None,
    seed: int = 1234,
    n_jobs: int = 1,
) -> PLSResult:
    """Run PLS between an imaging map and atlas gene expression.

    Parameters
    ----------
    data
        Imaging input in the same formats accepted by :func:`run_corr`.
    atlas
        Atlas identifier registered in the packaged atlas registry.
    hemisphere
        Hemisphere subset to analyse.
    regions
        Atlas region subset. ``"default"`` and ``"cort"`` keep cortex only;
        ``"all"`` and ``"cort+sub"`` include the packaged aseg add-on.
    source_space
        Declared space of the input image when resampling is required.
    input_rh
        Optional right-hemisphere surface file when ``data`` points to the left
        hemisphere file.
    n_components
        Number of PLS components to retain. Supply this or ``var``.
    var
        Alternative stopping rule based on cumulative explained variance.
    n_permutations
        Number of spatial null permutations used for component and gene-wise
        inference.
    null_method
        Spatial null model to use. ``"auto"`` chooses a supported method based
        on the atlas and available dependencies.
    output_dir
        Optional output directory for TSV, metadata, and plot files.
    enrichment_method
        Enrichment backend to apply to each retained component. Use
        ``"ensemble"``, ``"gsea"``, ``"ora"``, or ``"none"``. When omitted,
        the default is ``"ensemble"`` unless the legacy ``run_gsea`` or
        ``ora_p_threshold`` arguments imply another choice.
    run_gsea
        Legacy compatibility flag that forces the preranked GSEA backend.
    gene_set
        Packaged geneset name, Enrichr library name, or GMT path used by GSEA
        and ORA.
    geneset_organism
        Organism used when ``gene_set`` refers to an Enrichr/GSEApy library
        name rather than a packaged or local GMT file.
    ora_p_threshold
        If provided, run ORA on genes whose raw p-value is at or below this
        threshold, split into positive and negative component loadings.
    seed
        Random seed for the spatial null generation.
    n_jobs
        Number of worker processes used for PLS permutation fitting.

    Returns
    -------
    PLSResult
        Structured result object containing metadata, regional values, per-
        component statistics, and optional enrichment outputs.
    """

    config = build_run_config(
        "pls",
        atlas=atlas,
        hemisphere=hemisphere,
        regions=regions,
        source_space=source_space,
        n_permutations=n_permutations,
        null_method=null_method,
        output_dir=output_dir,
        enrichment_method=enrichment_method,
        run_gsea=run_gsea,
        gene_set=gene_set,
        geneset_organism=geneset_organism,
        ora_p_threshold=ora_p_threshold,
        n_components=n_components,
        var=var,
        seed=seed,
        n_jobs=n_jobs,
    )
    return _run_pls_configured(data, config, input_rh=input_rh)


__all__ = [
    "RunConfig",
    "build_run_config",
    "run_gene",
    "run_analysis",
    "run_corr",
    "run_pls",
]


def run_gene(
    gene: str,
    *,
    atlas: str = "dk",
    hemisphere: str = "left",
    regions: str = "default",
    zscore_expression: bool = True,
    top_n: int = 25,
    fdr_threshold: float = 0.05,
    output_dir=None,
):
    """Return one gene's regional expression vector and top co-expressed genes."""

    return _run_gene(
        gene,
        atlas=atlas,
        hemisphere=hemisphere,
        regions=regions,
        zscore_expression=zscore_expression,
        top_n=top_n,
        fdr_threshold=fdr_threshold,
        output_dir=output_dir,
    )
