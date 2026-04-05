"""Workflow implementation for GEDAR-style weighted gene-expression scoring."""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from ..exceptions import ConfigurationError
from ..gsea_utils import make_prerank_table, result_column, result_terms, run_prerank
from ..models import (
    AtlasSelection,
    ExpressionNormalization,
    GEDARDirection,
    GEDARResult,
    HemisphereMode,
    RankMode,
    RegionScope,
    WeightNormalization,
)
from ..genesets import resolve_geneset_resource
from ..ora import ora_from_gene_table
from ..preprocessing.gene_weights import (
    apply_brain_gene_filter_to_weights,
    apply_rank_selection_mask,
    clean_weight_table,
    load_weights_table,
    match_weight_table_to_genes,
    validate_rank_selection_args,
)
from ..scoring.weighted_expression import direction_mask, direction_modes, project_weighted_expression


def _resolve_gedar_enrichment_method(
    enrichment_method: str | None,
    *,
    run_gsea: bool,
) -> str:
    """Resolve the requested GEDAR enrichment mode.

    GEDAR keeps enrichment disabled by default. Users must explicitly request
    ``gsea`` or ``ora``.
    """

    if enrichment_method is not None:
        if enrichment_method not in {"gsea", "ora", "none"}:
            raise ConfigurationError(
                "GEDAR only supports enrichment_method='gsea', 'ora', or 'none'."
            )
        return enrichment_method
    if run_gsea:
        return "gsea"
    return "none"


def _run_gedar_gsea(
    gene_table: pd.DataFrame,
    *,
    gene_set: str,
    geneset_organism: str,
    gene_limit: int = 1500,
    n_perm: int = 1_000,
) -> pd.DataFrame:
    """Run preranked GSEA on the full matched signed GEDAR gene table."""

    try:
        import gseapy
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError("gseapy is required to run GSEA analyses.") from exc

    resource = resolve_geneset_resource(gene_set, organism=geneset_organism)
    ranking = make_prerank_table(
        gene_table["gene"].astype(str).tolist(),
        gene_table["input_weight"].to_numpy(dtype=float),
    )
    results = run_prerank(
        gseapy,
        ranking,
        resource,
        max_size=gene_limit,
        outdir=None,
        permutation_num=n_perm,
        seed=1234,
    )
    return pd.DataFrame.from_dict(
        OrderedDict(
            Term=result_terms(results.res2d),
            es=result_column(results.res2d, "ES", "es").astype(float),
            nes=result_column(results.res2d, "NES", "nes").astype(float),
            p_val=result_column(results.res2d, "NOM p-val", "pval", "p_val").astype(float),
            fdr=result_column(results.res2d, "FDR q-val", "fdr").astype(float),
        )
    )


def _run_gedar_ora(
    gene_table: pd.DataFrame,
    *,
    direction: GEDARDirection,
    gene_set: str,
    geneset_organism: str,
) -> dict[str, pd.DataFrame]:
    """Run ORA on the genes selected by GEDAR, split into up and down sets."""

    if direction == "split":
        selected = (
            gene_table["selected_up"].astype(bool).to_numpy(copy=False)
            | gene_table["selected_down"].astype(bool).to_numpy(copy=False)
        )
    else:
        selected = gene_table["selected"].astype(bool).to_numpy(copy=False)

    ora_table = gene_table.loc[:, ["gene", "input_weight"]].copy()
    ora_table["p_value"] = np.where(selected, 0.0, 1.0)
    return ora_from_gene_table(
        ora_table,
        gene_set=gene_set,
        geneset_organism=geneset_organism,
        score_column="input_weight",
        p_threshold=0.5,
    )


def run_gedar(
    weights,
    *,
    atlas: str = "dk",
    hemisphere: HemisphereMode = "both",
    regions: RegionScope = "default",
    gene_column: str = "gene",
    weight_column: str = "weight",
    rank_column: str | None = None,
    rank_mode: RankMode = "ascending",
    top_percent: float | None = None,
    top_n: int | None = None,
    p_threshold: float | None = None,
    direction: GEDARDirection = "combined",
    normalize_expression: ExpressionNormalization = "zscore",
    normalize_weights: WeightNormalization = "none",
    enrichment_method: str | None = None,
    run_gsea: bool = False,
    gene_set: str = "lake",
    geneset_organism: str = "Human",
    output_dir: str | Path | None = None,
    select_atlas_data_fn: Callable[..., AtlasSelection],
    load_brain_gene_symbols_fn: Callable[[], tuple[str, ...]],
    write_result_bundle_fn: Callable[[GEDARResult, Path], None] | None = None,
) -> GEDARResult:
    """Compute a PTRS-style weighted gene-expression score on a packaged atlas."""

    validate_rank_selection_args(
        rank_column=rank_column,
        top_percent=top_percent,
        top_n=top_n,
        p_threshold=p_threshold,
    )
    if normalize_expression not in {"zscore", "none"}:
        raise ConfigurationError("normalize_expression must be either 'zscore' or 'none'.")
    resolved_enrichment_method = _resolve_gedar_enrichment_method(
        enrichment_method,
        run_gsea=run_gsea,
    )

    weight_table, weights_source = load_weights_table(weights)
    if gene_column not in weight_table.columns:
        raise ConfigurationError(f"Column '{gene_column}' is not present in the weight table.")
    if weight_column not in weight_table.columns:
        raise ConfigurationError(f"Column '{weight_column}' is not present in the weight table.")
    if rank_column is not None and rank_column not in weight_table.columns:
        raise ConfigurationError(f"Column '{rank_column}' is not present in the weight table.")

    cleaned = clean_weight_table(
        weight_table,
        gene_column=gene_column,
        weight_column=weight_column,
        rank_column=rank_column,
        rank_mode=rank_mode,
    )
    if cleaned.table.empty:
        raise ConfigurationError(
            "No valid genes remained after filtering missing, invalid, and duplicated rows in the weight table."
        )

    filtered = apply_brain_gene_filter_to_weights(
        cleaned,
        gene_column=gene_column,
        brain_gene_symbols=set(load_brain_gene_symbols_fn()),
    )
    if filtered.table.empty:
        raise ConfigurationError("No genes remained after applying the packaged AHPA brain-gene filter.")

    selection = select_atlas_data_fn(
        atlas=atlas,
        hemisphere=hemisphere,
        regions=regions,
        zscore_expression=(normalize_expression == "zscore"),
    )
    available_genes = selection.expression.columns[2:].to_numpy(dtype=str, copy=False)
    matched = match_weight_table_to_genes(
        filtered.table,
        gene_column=gene_column,
        weight_column=weight_column,
        rank_column=rank_column,
        available_genes=available_genes,
    )
    if matched.table.empty:
        raise ConfigurationError(
            f"None of the requested genes were found in atlas '{selection.atlas.id}'."
        )

    expression_matrix = selection.expression.loc[:, matched.matched_genes].to_numpy(dtype=float, copy=True)
    regional_payload: dict[str, np.ndarray] = {}
    any_selected = False

    input_weights = matched.table["input_weight"].to_numpy(dtype=float)
    direction_columns: list[str] = []
    for mode in direction_modes(direction):
        mask = direction_mask(input_weights, mode)
        selected = np.zeros(matched.table.shape[0], dtype=bool)
        selected[np.flatnonzero(mask)] = apply_rank_selection_mask(
            matched.table.loc[mask].reset_index(drop=True),
            rank_mode=rank_mode,
            top_percent=top_percent,
            top_n=top_n,
            p_threshold=p_threshold,
        )
        projection = project_weighted_expression(
            expression_matrix,
            input_weights,
            selected,
            direction=mode,
            normalize_weights=normalize_weights,
            n_regions=selection.n_regions,
        )
        any_selected = any_selected or bool(np.any(projection.selected))

        if direction == "split":
            matched.table[f"weight_{mode}"] = projection.weights
            matched.table[f"selected_{mode}"] = projection.selected
            regional_payload[f"score_{mode}"] = projection.score
            regional_payload[f"score_{mode}_z"] = projection.score_z
            direction_columns.extend([f"weight_{mode}", f"selected_{mode}"])
        else:
            matched.table["weight"] = projection.weights
            matched.table["selected"] = projection.selected
            regional_payload["score"] = projection.score
            regional_payload["score_z"] = projection.score_z
            direction_columns.extend(["weight", "selected"])

    if not any_selected:
        raise ConfigurationError("No genes remained after applying the requested GEDAR filters.")

    regional_scores = selection.labels.reset_index(drop=True).assign(**regional_payload)
    matched.table["direction"] = direction
    ordered_columns = ["gene", "input_gene", "input_weight"]
    if "rank_value" in matched.table.columns:
        ordered_columns.append("rank_value")
    ordered_columns.extend(direction_columns)
    ordered_columns.append("direction")
    gene_table = matched.table.loc[:, ordered_columns].copy()
    gsea_table = None
    ora_tables = None
    if resolved_enrichment_method == "gsea":
        gsea_table = _run_gedar_gsea(
            gene_table,
            gene_set=gene_set,
            geneset_organism=geneset_organism,
        )
    elif resolved_enrichment_method == "ora":
        ora_tables = _run_gedar_ora(
            gene_table,
            direction=direction,
            gene_set=gene_set,
            geneset_organism=geneset_organism,
        )

    result = GEDARResult(
        atlas_id=selection.atlas.id,
        atlas_label=selection.atlas.label,
        hemisphere=selection.hemisphere,
        regions=selection.regions,
        requested_genes=filtered.requested_genes,
        regional_scores=regional_scores,
        gene_table=gene_table,
        excluded_table=filtered.excluded_table,
        gsea_table=gsea_table,
        ora_tables=ora_tables,
        matched_genes=matched.matched_genes,
        missing_genes=matched.missing_genes,
        weights_source=weights_source,
        gene_column=gene_column,
        weight_column=weight_column,
        rank_column=rank_column,
        rank_mode=rank_mode,
        direction=direction,
        normalize_expression=normalize_expression,
        normalize_weights=normalize_weights,
        top_percent=top_percent,
        top_n=top_n,
        p_threshold=p_threshold,
        enrichment_method=resolved_enrichment_method,
        geneset=gene_set if resolved_enrichment_method != "none" else None,
        geneset_organism=geneset_organism if resolved_enrichment_method != "none" else None,
        output_dir=None if output_dir is None else Path(output_dir),
    )
    if output_dir is not None and write_result_bundle_fn is not None:
        write_result_bundle_fn(result, Path(output_dir))
    return result
