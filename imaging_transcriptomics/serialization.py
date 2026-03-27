from __future__ import annotations

import json
import tempfile
from dataclasses import asdict
from pathlib import Path
from typing import Callable
import warnings

import numpy as np
import pandas as pd

from .exceptions import PlottingUnavailableError
from .models import CorrelationResult, GenePCAResult, PLSResult
from .text import write_readme


def metadata_dict(result: CorrelationResult | PLSResult | GenePCAResult) -> dict[str, object]:
    if isinstance(result, GenePCAResult):
        return {
            "method": "gene-pca",
            "atlas_id": result.atlas_id,
            "atlas_label": result.atlas_label,
            "hemisphere": result.hemisphere,
            "regions": result.regions,
            "n_genes_requested": len(result.requested_genes),
            "n_genes_used": len(result.matched_genes),
            "n_genes_missing": len(result.missing_genes),
            "n_components": int(result.variance_table.shape[0]),
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


def write_metadata_json(result: CorrelationResult | PLSResult | GenePCAResult, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "metadata.json"
    path.write_text(json.dumps(metadata_dict(result), indent=2, sort_keys=True) + "\n")
    return path


def _write_optional_plots(result: CorrelationResult | PLSResult | GenePCAResult, output_dir: Path) -> bool:
    from .plotting import save_result_plots

    try:
        save_result_plots(result, output_dir)
    except PlottingUnavailableError as exc:
        warnings.warn(str(exc), RuntimeWarning, stacklevel=2)
        return False
    return True


def _write_analysis_bundle(result: CorrelationResult | PLSResult, output_dir: Path, *, plots_available: bool) -> None:
    result.regional_values.to_csv(output_dir / "regional_values.tsv", sep="\t", index=False)
    write_metadata_json(result, output_dir)
    write_readme(result, output_dir, plots_available=plots_available)
    if isinstance(result, CorrelationResult):
        result.gene_table.to_csv(output_dir / "corr_genes.tsv", sep="\t", index=False)
        if result.gsea_table is not None:
            result.gsea_table.to_csv(output_dir / "gsea_corr_results.tsv", sep="\t", index=False)
        if result.ora_tables is not None:
            for direction, table in result.ora_tables.items():
                if table is not None:
                    table.to_csv(output_dir / f"ora_corr_{direction}.tsv", sep="\t", index=False)
        return

    pd.DataFrame(
        {
            "component": [component.index for component in result.components],
            "variance_explained": [component.explained_variance for component in result.components],
            "cumulative_variance": np.asarray(result.cumulative_variance[: len(result.components)], dtype=float),
            "permutation_p_value": [component.p_value for component in result.components],
        }
    ).to_csv(output_dir / "pls_summary.tsv", sep="\t", index=False)
    for component in result.components:
        component.gene_table.to_csv(output_dir / f"pls_component_{component.index}.tsv", sep="\t", index=False)
        if component.gsea_table is not None:
            component.gsea_table.to_csv(
                output_dir / f"gsea_pls{component.index}_results.tsv",
                sep="\t",
                index=False,
            )
        if component.ora_tables is not None:
            for direction, table in component.ora_tables.items():
                if table is not None:
                    table.to_csv(
                        output_dir / f"ora_pls{component.index}_{direction}.tsv",
                        sep="\t",
                        index=False,
                    )


def _write_gene_pca_bundle(result: GenePCAResult, output_dir: Path, *, plots_available: bool) -> None:
    result.regional_scores.to_csv(output_dir / "gene_pca_scores.tsv", sep="\t", index=False)
    result.gene_loadings.to_csv(output_dir / "gene_pca_loadings.tsv", sep="\t", index=False)
    result.variance_table.to_csv(output_dir / "gene_pca_variance.tsv", sep="\t", index=False)
    (output_dir / "matched_genes.txt").write_text("\n".join(result.matched_genes) + ("\n" if result.matched_genes else ""))
    (output_dir / "missing_genes.txt").write_text("\n".join(result.missing_genes) + ("\n" if result.missing_genes else ""))
    write_metadata_json(result, output_dir)
    write_readme(result, output_dir, plots_available=plots_available)


def write_result_bundle(result: CorrelationResult | PLSResult | GenePCAResult, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_available = _write_optional_plots(result, output_dir)
    if isinstance(result, GenePCAResult):
        _write_gene_pca_bundle(result, output_dir, plots_available=plots_available)
        return
    _write_analysis_bundle(result, output_dir, plots_available=plots_available)


def read_optional_table(path: Path) -> pd.DataFrame | None:
    return pd.read_csv(path, sep="\t") if path.exists() else None


def run_to_tables(
    output_dir: Path | None,
    runner: Callable[[Path], None],
    table_paths: Callable[[Path], list[Path]],
) -> list[pd.DataFrame | None]:
    if output_dir is not None:
        runner(output_dir)
        return [read_optional_table(path) for path in table_paths(output_dir)]

    with tempfile.TemporaryDirectory() as tmpdir:
        temp_dir = Path(tmpdir)
        runner(temp_dir)
        return [read_optional_table(path) for path in table_paths(temp_dir)]
