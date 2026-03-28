"""Metadata helpers for persisted analysis result bundles."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from ..models import CorrelationResult, GEDARResult, GenePCAResult, PLSResult


def metadata_dict(result: CorrelationResult | PLSResult | GenePCAResult | GEDARResult) -> dict[str, object]:
    """Return the JSON-serializable metadata payload for one result object."""

    if isinstance(result, GenePCAResult):
        return {
            "method": "gene-pca",
            "atlas_id": result.atlas_id,
            "atlas_label": result.atlas_label,
            "hemisphere": result.hemisphere,
            "regions": result.regions,
            "brain_gene_filter": "AHPA_mrna_brain.tsv",
            "n_genes_requested": len(result.requested_genes),
            "n_genes_used": len(result.matched_genes),
            "n_genes_filtered": len(result.brain_filtered_genes),
            "n_genes_missing": len(result.missing_genes),
            "n_components": int(result.variance_table.shape[0]),
            "requested_genes": list(result.requested_genes),
            "matched_genes": list(result.matched_genes),
            "brain_filtered_genes": list(result.brain_filtered_genes),
            "missing_genes": list(result.missing_genes),
            "output_dir": None if result.output_dir is None else str(result.output_dir),
        }
    if isinstance(result, GEDARResult):
        exclusion_counts = (
            result.excluded_table["exclusion_reason"].astype(str).value_counts().sort_index().to_dict()
            if not result.excluded_table.empty and "exclusion_reason" in result.excluded_table.columns
            else {}
        )
        if result.direction == "split":
            n_selected = {
                "up": int(result.gene_table["selected_up"].astype(bool).sum()),
                "down": int(result.gene_table["selected_down"].astype(bool).sum()),
            }
        else:
            n_selected = int(result.gene_table["selected"].astype(bool).sum())
        return {
            "method": "gedar",
            "atlas_id": result.atlas_id,
            "atlas_label": result.atlas_label,
            "hemisphere": result.hemisphere,
            "regions": result.regions,
            "weights_source": result.weights_source,
            "gene_column": result.gene_column,
            "weight_column": result.weight_column,
            "rank_column": result.rank_column,
            "rank_mode": result.rank_mode,
            "direction": result.direction,
            "normalize_expression": result.normalize_expression,
            "normalize_weights": result.normalize_weights,
            "brain_gene_filter": "AHPA_mrna_brain.tsv",
            "top_percent": result.top_percent,
            "top_n": result.top_n,
            "p_threshold": result.p_threshold,
            "n_genes_requested": len(result.requested_genes),
            "n_genes_used": len(result.matched_genes),
            "n_genes_missing": len(result.missing_genes),
            "n_genes_selected": n_selected,
            "n_rows_excluded": int(result.excluded_table.shape[0]),
            "excluded_counts": exclusion_counts,
            "requested_genes": list(result.requested_genes),
            "matched_genes": list(result.matched_genes),
            "missing_genes": list(result.missing_genes),
            "output_dir": None if result.output_dir is None else str(result.output_dir),
        }

    metadata = asdict(result.metadata)
    metadata["output_dir"] = None if result.output_dir is None else str(result.output_dir)
    if isinstance(result, PLSResult):
        metadata["cumulative_variance"] = np.asarray(result.cumulative_variance, dtype=float).tolist()
    return metadata


def write_metadata_json(result: CorrelationResult | PLSResult | GenePCAResult | GEDARResult, output_dir: Path) -> Path:
    """Write the metadata payload for one result bundle as JSON."""

    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "metadata.json"
    path.write_text(json.dumps(metadata_dict(result), indent=2, sort_keys=True) + "\n")
    return path
