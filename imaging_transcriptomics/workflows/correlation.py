from __future__ import annotations

from ..corr import CorrAnalysis
from ..models import CorrelationResult
from ..nulls import permute_scan_values as _default_permute_scan_values
from ..scan import regional_values_frame
from ..serialization import run_to_tables, write_result_bundle
from .shared import corr_gene_table, prepare_analysis_inputs, result_metadata


def run_corr_configured(
    data,
    config,
    *,
    input_rh=None,
    permute_scan_values_fn=_default_permute_scan_values,
) -> CorrelationResult:
    """Execute a configured spatial correlation workflow."""

    extracted, gene_exp, gene_labels, imaging = prepare_analysis_inputs(data, config, input_rh=input_rh)
    permuted, resolved_null_method = permute_scan_values_fn(
        extracted,
        n_permutations=config.n_permutations,
        null_method=config.null_method,
        seed=config.seed,
    )

    analysis = CorrAnalysis(
        n_iterations=config.n_permutations,
        n_genes=gene_labels.shape[0],
        store_boot_corr=config.enrichment_method in {"gsea", "ensemble"},
        n_jobs=config.n_jobs,
    )
    analysis.bootstrap_correlation(imaging, permuted, gene_exp, gene_labels)

    gsea_table = None
    ensemble_table = None
    ora_tables = None
    if config.enrichment_method == "gsea":
        gsea_table = run_to_tables(
            config.output_dir,
            lambda outdir: analysis.gsea(
                gene_set=config.gene_set,
                outdir=outdir,
                geneset_organism=config.geneset_organism,
            ),
            lambda outdir: [outdir / "gsea_corr_results.tsv"],
        )[0]
    elif config.enrichment_method == "ensemble":
        ensemble_table = run_to_tables(
            config.output_dir,
            lambda outdir: analysis.ensemble(
                gene_set=config.gene_set,
                outdir=outdir,
                geneset_organism=config.geneset_organism,
            ),
            lambda outdir: [outdir / "ensemble_corr_results.tsv"],
        )[0]
    elif config.enrichment_method == "ora":
        ora_up, ora_down = run_to_tables(
            config.output_dir,
            lambda outdir: analysis.ora(
                gene_set=config.gene_set,
                outdir=outdir,
                p_threshold=config.ora_p_threshold,
                geneset_organism=config.geneset_organism,
            ),
            lambda outdir: [outdir / "ora_corr_up.tsv", outdir / "ora_corr_down.tsv"],
        )
        ora_tables = {"up": ora_up, "down": ora_down}

    result = CorrelationResult(
        metadata=result_metadata(extracted, config, null_method=resolved_null_method),
        regional_values=regional_values_frame(extracted),
        gene_table=corr_gene_table(analysis),
        gsea_table=gsea_table,
        ensemble_table=ensemble_table,
        ora_tables=ora_tables,
        output_dir=config.output_dir,
    )
    if config.output_dir is not None:
        write_result_bundle(result, config.output_dir)
    return result
