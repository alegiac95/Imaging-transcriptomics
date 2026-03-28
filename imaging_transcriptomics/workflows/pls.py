from __future__ import annotations

import numpy as np

from ..models import PLSResult
from ..nulls import permute_scan_values as _default_permute_scan_values
from ..pls import PLSAnalysis
from ..scan import regional_values_frame
from ..serialization import run_to_tables, write_result_bundle
from .shared import pls_components, prepare_analysis_inputs, result_metadata


def run_pls_configured(
    data,
    config,
    *,
    input_rh=None,
    permute_scan_values_fn=_default_permute_scan_values,
) -> PLSResult:
    """Execute a configured PLS workflow."""

    extracted, gene_exp, gene_labels, imaging = prepare_analysis_inputs(data, config, input_rh=input_rh)
    permuted, resolved_null_method = permute_scan_values_fn(
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
    analysis.gene_results.results.compute(n_jobs=config.n_jobs)

    if config.run_gsea:
        gsea_tables = run_to_tables(
            config.output_dir,
            lambda outdir: analysis.gene_results.results.gsea(
                gene_set=config.gene_set,
                outdir=outdir,
                geneset_organism=config.geneset_organism,
            ),
            lambda outdir: [
                outdir / f"gsea_pls{index}_results.tsv"
                for index in range(1, analysis.n_components + 1)
            ],
        )
    else:
        gsea_tables = [None] * analysis.n_components

    if config.ora_p_threshold is not None:
        ora_files = run_to_tables(
            config.output_dir,
            lambda outdir: analysis.gene_results.results.ora(
                gene_set=config.gene_set,
                outdir=outdir,
                p_threshold=config.ora_p_threshold,
                geneset_organism=config.geneset_organism,
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
        metadata=result_metadata(
            extracted,
            config,
            null_method=resolved_null_method,
            n_components=analysis.n_components,
        ),
        regional_values=regional_values_frame(extracted),
        components=pls_components(analysis, gsea_tables, ora_tables),
        cumulative_variance=np.cumsum(analysis.components_var),
        output_dir=config.output_dir,
    )
    if config.output_dir is not None:
        write_result_bundle(result, config.output_dir)
    return result
