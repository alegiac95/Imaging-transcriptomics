"""Writers for persisted result bundles and temporary table extraction."""

from __future__ import annotations

import tempfile
import warnings
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from ..exceptions import PlottingUnavailableError
from ..models import CorrelationResult, GEDARResult, GenePCAResult, PLSResult
from ..text import write_readme
from .metadata import write_metadata_json


def _write_optional_plots(result: CorrelationResult | PLSResult | GenePCAResult | GEDARResult, output_dir: Path) -> bool:
    """Try to write plots for one result bundle and warn on optional failures."""

    from ..plotting import save_result_plots

    try:
        save_result_plots(result, output_dir)
    except PlottingUnavailableError as exc:
        warnings.warn(str(exc), RuntimeWarning, stacklevel=2)
        return False
    return True


def _write_analysis_bundle(result: CorrelationResult | PLSResult, output_dir: Path, *, plots_available: bool) -> None:
    """Write the standard persisted bundle for correlation and PLS workflows."""

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
            "p": [component.p_value for component in result.components],
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
    """Write the standard persisted bundle for gene-PCA runs."""

    result.regional_scores.to_csv(output_dir / "gene_pca_scores.tsv", sep="\t", index=False)
    result.gene_loadings.to_csv(output_dir / "gene_pca_loadings.tsv", sep="\t", index=False)
    result.variance_table.to_csv(output_dir / "gene_pca_variance.tsv", sep="\t", index=False)
    (output_dir / "matched_genes.txt").write_text("\n".join(result.matched_genes) + ("\n" if result.matched_genes else ""))
    (output_dir / "brain_filtered_genes.txt").write_text("\n".join(result.brain_filtered_genes) + ("\n" if result.brain_filtered_genes else ""))
    (output_dir / "missing_genes.txt").write_text("\n".join(result.missing_genes) + ("\n" if result.missing_genes else ""))
    write_metadata_json(result, output_dir)
    write_readme(result, output_dir, plots_available=plots_available)


def _write_gedar_bundle(result: GEDARResult, output_dir: Path, *, plots_available: bool) -> None:
    """Write the standard persisted bundle for GEDAR runs."""

    result.regional_scores.to_csv(output_dir / "gedar_scores.tsv", sep="\t", index=False)
    result.gene_table.to_csv(output_dir / "gedar_genes.tsv", sep="\t", index=False)
    result.excluded_table.to_csv(output_dir / "gedar_excluded.tsv", sep="\t", index=False)
    (output_dir / "matched_genes.txt").write_text("\n".join(result.matched_genes) + ("\n" if result.matched_genes else ""))
    (output_dir / "missing_genes.txt").write_text("\n".join(result.missing_genes) + ("\n" if result.missing_genes else ""))
    write_metadata_json(result, output_dir)
    write_readme(result, output_dir, plots_available=plots_available)


def write_result_bundle(result: CorrelationResult | PLSResult | GenePCAResult | GEDARResult, output_dir: Path) -> None:
    """Write all standard tables, metadata, README text, and optional plots."""

    output_dir.mkdir(parents=True, exist_ok=True)
    plots_available = _write_optional_plots(result, output_dir)
    if isinstance(result, GenePCAResult):
        _write_gene_pca_bundle(result, output_dir, plots_available=plots_available)
        return
    if isinstance(result, GEDARResult):
        _write_gedar_bundle(result, output_dir, plots_available=plots_available)
        return
    _write_analysis_bundle(result, output_dir, plots_available=plots_available)


def read_optional_table(path: Path) -> pd.DataFrame | None:
    """Read one TSV table if it exists, otherwise return ``None``."""

    return pd.read_csv(path, sep="\t") if path.exists() else None


def run_to_tables(
    output_dir: Path | None,
    runner: Callable[[Path], None],
    table_paths: Callable[[Path], list[Path]],
) -> list[pd.DataFrame | None]:
    """Run one side-effecting writer and read back the requested tables."""

    if output_dir is not None:
        runner(output_dir)
        return [read_optional_table(path) for path in table_paths(output_dir)]

    with tempfile.TemporaryDirectory() as tmpdir:
        temp_dir = Path(tmpdir)
        runner(temp_dir)
        return [read_optional_table(path) for path in table_paths(temp_dir)]
