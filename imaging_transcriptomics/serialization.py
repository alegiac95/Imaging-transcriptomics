from __future__ import annotations

import json
import tempfile
from dataclasses import asdict
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .models import CorrelationResult, PLSResult
from .text import write_readme


def metadata_dict(result: CorrelationResult | PLSResult) -> dict[str, object]:
    metadata = asdict(result.metadata)
    metadata["output_dir"] = None if result.output_dir is None else str(result.output_dir)
    if isinstance(result, PLSResult):
        metadata["cumulative_variance"] = np.asarray(result.cumulative_variance, dtype=float).tolist()
    return metadata


def write_metadata_json(result: CorrelationResult | PLSResult, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "metadata.json"
    path.write_text(json.dumps(metadata_dict(result), indent=2, sort_keys=True) + "\n")
    return path


def write_result_bundle(result: CorrelationResult | PLSResult, output_dir: Path) -> None:
    from .plotting import save_result_plots

    output_dir.mkdir(parents=True, exist_ok=True)
    result.regional_values.to_csv(output_dir / "regional_values.tsv", sep="\t", index=False)
    write_metadata_json(result, output_dir)
    save_result_plots(result, output_dir)
    write_readme(result, output_dir)
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
