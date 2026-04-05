from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import RunConfig
from ..gene_expression import expression_matrix
from ..models import AnalysisMetadata, PLSComponentResult
from ..scan import extract_scan_data


def standardize_vector(values: np.ndarray) -> np.ndarray:
    """Center and scale a one-dimensional vector with sample standard deviation."""

    vector = np.asarray(values, dtype=float).reshape(-1)
    centered = vector - vector.mean()
    scale = vector.std(ddof=1)
    if scale == 0:
        return centered
    return centered / scale


def prepare_analysis_inputs(data, config: RunConfig, *, input_rh=None):
    """Extract atlas-aligned imaging and expression inputs for one run config."""

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
        standardize_vector(extracted.values),
    )


def result_metadata(extracted, config: RunConfig, *, null_method: str, n_components: int | None = None) -> AnalysisMetadata:
    """Build the shared metadata payload for correlation and PLS outputs."""

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
        enrichment_method=config.enrichment_method,
        geneset=config.gene_set if config.enrichment_method != "none" else None,
        geneset_organism=config.geneset_organism if config.enrichment_method != "none" else None,
        ora_p_threshold=config.ora_p_threshold if config.enrichment_method == "ora" else None,
        n_components=n_components,
    )


def corr_gene_table(analysis) -> pd.DataFrame:
    """Convert correlation gene statistics into the public result schema."""

    return pd.DataFrame(
        {
            "gene": analysis.gene_results.results.genes[:, 0],
            "score": analysis.gene_results.results.corr[0, :],
            "p": analysis.gene_results.results.pval[0, :],
            "fdr": analysis.gene_results.results.pval_corr[0, :],
            "maxT": analysis.gene_results.results.pval_fwer[0, :],
        }
    )


def pls_components(
    analysis,
    gsea_tables: list[pd.DataFrame | None],
    ensemble_tables: list[pd.DataFrame | None],
    ora_tables: list[dict[str, pd.DataFrame] | None],
) -> tuple[PLSComponentResult, ...]:
    """Pack per-component PLS outputs into typed result records."""

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
                    "p": analysis.gene_results.results.boot.pval[index, :],
                    "fdr": analysis.gene_results.results.boot.pval_corr[index, :],
                    "maxT": analysis.gene_results.results.boot.pval_fwer[index, :],
                }
            ),
            gsea_table=gsea_tables[index],
            ensemble_table=ensemble_tables[index],
            ora_tables=ora_tables[index],
        )
        for index in range(analysis.n_components)
    )
