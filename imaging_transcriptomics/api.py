from __future__ import annotations

import numpy as np
import pandas as pd

from .config import (
    DEFAULT_NULL_METHOD,
    DEFAULT_PERMUTATIONS,
    RunConfig,
    build_run_config,
)
from .gene_expression import expression_matrix
from .models import AnalysisMetadata, CorrelationResult, PLSComponentResult, PLSResult
from .nulls import (
    permute_scan_values as _permute_scan_values,
)
from .scan import extract_scan_data, regional_values_frame
from .serialization import run_to_tables as _run_to_tables, write_result_bundle


def _standardize_vector(values: np.ndarray) -> np.ndarray:
    vector = np.asarray(values, dtype=float).reshape(-1)
    centered = vector - vector.mean()
    scale = vector.std(ddof=1)
    if scale == 0:
        return centered
    return centered / scale


def _prepare_analysis_inputs(data, config: RunConfig, *, input_rh=None):
    extracted = extract_scan_data(
        data,
        atlas=config.atlas,
        hemisphere=config.hemisphere,
        regions=config.regions,
        source_space=config.source_space,
        input_rh=input_rh,
    )
    return (
        extracted,
        expression_matrix(extracted.selection),
        extracted.selection.gene_labels,
        _standardize_vector(extracted.values),
    )


def _metadata(extracted, config: RunConfig, *, null_method: str, n_components: int | None = None) -> AnalysisMetadata:
    return AnalysisMetadata(
        method=config.method,
        atlas_id=extracted.selection.atlas.id,
        atlas_label=extracted.selection.atlas.label,
        hemisphere=extracted.selection.hemisphere,
        regions=extracted.selection.regions,
        source=extracted.source,
        source_kind=extracted.source_kind,
        source_space=extracted.source_space,
        n_permutations=config.n_permutations,
        null_method=null_method,
        geneset=config.gene_set if (config.run_gsea or config.ora_p_threshold is not None) else None,
        ora_p_threshold=config.ora_p_threshold,
        n_components=n_components,
    )


def _corr_gene_table(analysis) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene": analysis.gene_results.results.genes[:, 0],
            "score": analysis.gene_results.results.corr[0, :],
            "p_value": analysis.gene_results.results.pval[0, :],
            "fdr": analysis.gene_results.results.pval_corr[0, :],
            "fwer_maxT": analysis.gene_results.results.pval_fwer[0, :],
        }
    )


def _pls_components(
    analysis,
    gsea_tables: list[pd.DataFrame | None],
    ora_tables: list[dict[str, pd.DataFrame] | None],
) -> tuple[PLSComponentResult, ...]:
    return tuple(
        PLSComponentResult(
            index=index + 1,
            explained_variance=float(analysis.components_var[index]),
            p_value=float(analysis.p_val[index]),
            gene_table=pd.DataFrame(
                {
                    "gene": analysis.gene_results.results.boot.genes[index, :],
                    "weight": analysis.gene_results.results.boot.weights_sorted[index, :],
                    "zscore": analysis.gene_results.results.boot.z_score[index, :],
                    "p_value": analysis.gene_results.results.boot.pval[index, :],
                    "fdr": analysis.gene_results.results.boot.pval_corr[index, :],
                    "fwer_maxT": analysis.gene_results.results.boot.pval_fwer[index, :],
                }
            ),
            gsea_table=gsea_tables[index],
            ora_tables=ora_tables[index],
        )
        for index in range(analysis.n_components)
    )


def run_analysis(data, config: RunConfig, *, input_rh=None) -> CorrelationResult | PLSResult:
    if config.method == "corr":
        return _run_corr_configured(data, config, input_rh=input_rh)
    return _run_pls_configured(data, config, input_rh=input_rh)


def _run_corr_configured(data, config: RunConfig, *, input_rh=None) -> CorrelationResult:
    from .corr import CorrAnalysis

    extracted, gene_exp, gene_labels, imaging = _prepare_analysis_inputs(data, config, input_rh=input_rh)
    permuted, resolved_null_method = _permute_scan_values(
        extracted,
        n_permutations=config.n_permutations,
        null_method=config.null_method,
        seed=config.seed,
    )

    analysis = CorrAnalysis(n_iterations=config.n_permutations, n_genes=gene_labels.shape[0])
    analysis.bootstrap_correlation(imaging, permuted, gene_exp, gene_labels)

    gsea_table = None
    if config.run_gsea:
        gsea_table = _run_to_tables(
            config.output_dir,
            lambda outdir: analysis.gsea(gene_set=config.gene_set, outdir=outdir),
            lambda outdir: [outdir / "gsea_corr_results.tsv"],
        )[0]

    ora_tables = None
    if config.ora_p_threshold is not None:
        ora_up, ora_down = _run_to_tables(
            config.output_dir,
            lambda outdir: analysis.ora(
                gene_set=config.gene_set,
                outdir=outdir,
                p_threshold=config.ora_p_threshold,
            ),
            lambda outdir: [outdir / "ora_corr_up.tsv", outdir / "ora_corr_down.tsv"],
        )
        ora_tables = {"up": ora_up, "down": ora_down}

    result = CorrelationResult(
        metadata=_metadata(extracted, config, null_method=resolved_null_method),
        regional_values=regional_values_frame(extracted),
        gene_table=_corr_gene_table(analysis),
        gsea_table=gsea_table,
        ora_tables=ora_tables,
        output_dir=config.output_dir,
    )
    if config.output_dir is not None:
        write_result_bundle(result, config.output_dir)
    return result


def _run_pls_configured(data, config: RunConfig, *, input_rh=None) -> PLSResult:
    from .pls import PLSAnalysis

    extracted, gene_exp, gene_labels, imaging = _prepare_analysis_inputs(data, config, input_rh=input_rh)
    permuted, resolved_null_method = _permute_scan_values(
        extracted,
        n_permutations=config.n_permutations,
        null_method=config.null_method,
        seed=config.seed,
    )

    analysis = PLSAnalysis(
        imaging,
        gene_exp,
        n_components=config.n_components,
        var=config.var,
        n_iter=config.n_permutations,
        n_jobs=config.n_jobs,
    )
    analysis.boot_pls(
        imaging,
        permuted,
        gene_exp,
        scan_data=extracted.values,
        gene_labels=gene_labels,
    )
    analysis.gene_results.results.compute()

    if config.run_gsea:
        gsea_tables = _run_to_tables(
            config.output_dir,
            lambda outdir: analysis.gene_results.results.gsea(gene_set=config.gene_set, outdir=outdir),
            lambda outdir: [
                outdir / f"gsea_pls{index}_results.tsv"
                for index in range(1, analysis.n_components + 1)
            ],
        )
    else:
        gsea_tables = [None] * analysis.n_components

    if config.ora_p_threshold is not None:
        ora_files = _run_to_tables(
            config.output_dir,
            lambda outdir: analysis.gene_results.results.ora(
                gene_set=config.gene_set,
                outdir=outdir,
                p_threshold=config.ora_p_threshold,
            ),
            lambda outdir: [
                outdir / f"ora_pls{index}_{direction}.tsv"
                for index in range(1, analysis.n_components + 1)
                for direction in ("up", "down")
            ],
        )
        ora_tables = [
            {"up": ora_files[2 * index], "down": ora_files[2 * index + 1]}
            for index in range(analysis.n_components)
        ]
    else:
        ora_tables = [None] * analysis.n_components

    result = PLSResult(
        metadata=_metadata(
            extracted,
            config,
            null_method=resolved_null_method,
            n_components=analysis.n_components,
        ),
        regional_values=regional_values_frame(extracted),
        components=_pls_components(analysis, gsea_tables, ora_tables),
        cumulative_variance=np.cumsum(analysis.components_var),
        output_dir=config.output_dir,
    )
    if config.output_dir is not None:
        write_result_bundle(result, config.output_dir)
    return result


def run_corr(
    data,
    *,
    atlas: str = "dk",
    hemisphere: str = "left",
    regions: str = "all",
    source_space: str | None = None,
    input_rh=None,
    n_permutations: int = DEFAULT_PERMUTATIONS,
    null_method: str = DEFAULT_NULL_METHOD,
    output_dir=None,
    run_gsea: bool = False,
    gene_set: str = "lake",
    ora_p_threshold: float | None = None,
    seed: int = 1234,
    n_jobs: int = 1,
) -> CorrelationResult:
    config = build_run_config(
        "corr",
        atlas=atlas,
        hemisphere=hemisphere,
        regions=regions,
        source_space=source_space,
        n_permutations=n_permutations,
        null_method=null_method,
        output_dir=output_dir,
        run_gsea=run_gsea,
        gene_set=gene_set,
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
    regions: str = "all",
    source_space: str | None = None,
    input_rh=None,
    n_components: int | None = None,
    var: float | None = None,
    n_permutations: int = DEFAULT_PERMUTATIONS,
    null_method: str = DEFAULT_NULL_METHOD,
    output_dir=None,
    run_gsea: bool = False,
    gene_set: str = "lake",
    ora_p_threshold: float | None = None,
    seed: int = 1234,
    n_jobs: int = 1,
) -> PLSResult:
    config = build_run_config(
        "pls",
        atlas=atlas,
        hemisphere=hemisphere,
        regions=regions,
        source_space=source_space,
        n_permutations=n_permutations,
        null_method=null_method,
        output_dir=output_dir,
        run_gsea=run_gsea,
        gene_set=gene_set,
        ora_p_threshold=ora_p_threshold,
        n_components=n_components,
        var=var,
        seed=seed,
        n_jobs=n_jobs,
    )
    return _run_pls_configured(data, config, input_rh=input_rh)
