from __future__ import annotations

import json
import inspect
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from .atlas_registry import get_atlas


MirrorMode = Literal[None, "bidirectional", "leftright", "rightleft"]


def _clean_labels(atlas_info) -> pd.DataFrame:
    labels = pd.read_csv(atlas_info)
    if "Unnamed: 0" in labels.columns:
        labels = labels.drop(columns=["Unnamed: 0"])
    return labels


def _expression_frame(expression: pd.DataFrame, atlas_info) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = expression.reset_index().rename(columns={expression.index.name or "index": "id"})
    labels = _clean_labels(atlas_info)
    gene_columns = [column for column in expression.columns if column not in {"id", "Region"}]
    aligned = labels[["id"]].merge(expression[["id", *gene_columns]], on="id", how="left", validate="1:1")
    if aligned[gene_columns].isna().any(axis=None):
        missing_ids = aligned.loc[aligned[gene_columns].isna().any(axis=1), "id"].astype(str).tolist()
        raise ValueError(f"Missing expression rows for atlas ids: {', '.join(missing_ids[:5])}")
    return labels, aligned


def _write_gene_labels(genes: np.ndarray, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    genes = np.asarray(genes, dtype=str)
    if path.exists():
        existing = np.load(path, allow_pickle=False).astype(str, copy=False)
        if existing.shape != genes.shape or np.any(existing != genes):
            raise ValueError(f"Shared gene labels file {path} does not match the current atlas gene ordering.")
        return
    np.save(path, genes)


def _write_expression_archive(values: np.ndarray, path: Path) -> None:
    np.savez_compressed(path, values=np.asarray(values, dtype=np.float32, copy=True))


def _patch_abagen_pandas_compat(abagen) -> None:
    probes = abagen.probes_
    io = abagen.io
    if getattr(probes._groupby_structure_id, "__name__", "") == "_groupby_structure_id_compat":
        patched_groupby = True
    else:
        patched_groupby = False

    if not patched_groupby:
        def _groupby_structure_id_compat(microarray, annotation):
            sid = io.read_annotation(annotation)["structure_id"]
            return io.read_microarray(microarray).T.groupby(sid).mean().T

        probes._groupby_structure_id = _groupby_structure_id_compat

    if "inplace" not in inspect.signature(pd.DataFrame.set_axis).parameters:
        original_set_axis = pd.DataFrame.set_axis
        if getattr(original_set_axis, "__name__", "") != "_set_axis_compat":
            def _set_axis_compat(self, labels, axis=0, inplace=None, copy=None):
                del inplace
                return original_set_axis(self, labels, axis=axis, copy=copy)

            pd.DataFrame.set_axis = _set_axis_compat


def _import_abagen():
    try:
        import abagen
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "abagen is required to build atlas expression assets. Install imaging-transcriptomics[maps]."
        ) from exc
    _patch_abagen_pandas_compat(abagen)
    return abagen


def build_expression_assets(
    atlas: str,
    output_dir,
    atlas_image=None,
    atlas_info=None,
    *,
    lr_mirror: MirrorMode = "leftright",
    missing: str | None = None,
    donors: str | list[str] = "all",
    n_proc: int = 1,
    data_dir=None,
    geometry=None,
    space: str | None = None,
) -> dict[str, Path]:
    abagen = _import_abagen()
    spec = get_atlas(atlas)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    atlas_image = atlas_image or spec.volume_1mm_path or spec.surface_paths
    atlas_info = atlas_info or spec.labels_path
    if isinstance(atlas_image, tuple):
        atlas_image = tuple(str(path) for path in atlas_image)
    if isinstance(geometry, tuple):
        geometry = tuple(str(path) for path in geometry)
    if atlas_image is None or atlas_info is None:
        raise ValueError(
            f"Atlas '{spec.id}' does not have enough packaged metadata to build expression assets automatically. "
            "Provide atlas_image and atlas_info explicitly."
        )
    if geometry is None and spec.surface_geometry is not None:
        geometry = spec.surface_geometry
    if space is None and spec.surface_space is not None and isinstance(atlas_image, tuple):
        space = spec.surface_space
    if geometry is not None:
        atlas_image = abagen.check_atlas(
            atlas_image,
            atlas_info=atlas_info,
            geometry=geometry,
            space=space,
            data_dir=data_dir,
        )

    expression, counts, report = abagen.get_expression_data(
        atlas_image,
        atlas_info,
        lr_mirror=lr_mirror,
        missing=missing,
        donors=donors,
        data_dir=data_dir,
        return_counts=True,
        return_report=True,
        n_proc=n_proc,
        verbose=1,
    )

    expression_name = (
        spec.expression_path.name if spec.expression_path is not None else f"atlas-{spec.id}_gene_expression_data.npz"
    )
    labels_name = spec.labels_path.name if spec.labels_path is not None else f"atlas-{spec.id}_labels.csv"
    expression_out = output_path / expression_name
    counts_out = output_path / f"atlas-{spec.id}_sample_counts.csv"
    report_out = output_path / "README.txt"
    provenance_out = output_path / "provenance.json"
    atlas_info_out = output_path / labels_name
    gene_labels_out = (
        spec.gene_labels_path
        if spec.gene_labels_path is not None
        else output_path / f"atlas-{spec.id}_gene_labels.npy"
    )

    labels, expression = _expression_frame(expression, atlas_info)
    genes = expression.columns[1:].to_numpy(dtype=str)
    _write_gene_labels(genes, gene_labels_out)
    _write_expression_archive(expression.iloc[:, 1:].to_numpy(dtype=np.float32, copy=True), expression_out)
    counts.to_csv(counts_out, index=True)
    if Path(atlas_info).exists():
        labels.to_csv(atlas_info_out, index=False)
    report_out.write_text(report)
    provenance_out.write_text(
        json.dumps(
            {
                "atlas": spec.id,
                "atlas_label": spec.label,
                "lr_mirror": lr_mirror,
                "missing": missing,
                "donors": donors,
                "n_proc": n_proc,
                "source_image": str(atlas_image),
                "source_info": str(atlas_info),
                "expression_format": "npz",
                "value_dtype": "float32",
                "gene_labels": str(gene_labels_out),
                "row_order": "aligned to labels csv order",
            },
            indent=2,
        )
    )
    return {
        "expression": expression_out,
        "gene_labels": gene_labels_out,
        "counts": counts_out,
        "report": report_out,
        "provenance": provenance_out,
        "labels": atlas_info_out,
    }
