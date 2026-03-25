from __future__ import annotations

import json
import inspect
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from .atlas_registry import get_atlas


MirrorMode = Literal[None, "bidirectional", "leftright", "rightleft"]


def _expression_frame(expression: pd.DataFrame, atlas_info) -> pd.DataFrame:
    expression = expression.reset_index().rename(columns={expression.index.name or "index": "id"})
    if "Region" in expression.columns:
        return expression

    labels = pd.read_csv(atlas_info)
    if "Unnamed: 0" in labels.columns:
        labels = labels.drop(columns=["Unnamed: 0"])
    regions = labels[["id", "label"]].copy()
    if "hemisphere" in labels.columns:
        hemi = labels["hemisphere"].astype(str)
        mask = hemi.isin(["L", "R"])
        regions["Region"] = labels["label"].astype(str)
        regions.loc[mask, "Region"] = hemi.loc[mask] + "_" + labels.loc[mask, "label"].astype(str)
    else:
        regions["Region"] = labels["label"].astype(str)
    return regions[["id", "Region"]].merge(expression, on="id", how="right")


def _write_expression_archive(expression: pd.DataFrame, path: Path) -> None:
    genes = expression.columns[2:].to_numpy(dtype=str)
    values = expression.iloc[:, 2:].to_numpy(dtype=np.float32, copy=True)
    np.savez_compressed(
        path,
        ids=expression.iloc[:, 0].to_numpy(dtype=np.int32, copy=True),
        regions=expression.iloc[:, 1].astype(str).to_numpy(dtype=str, copy=True),
        genes=genes,
        values=values,
    )


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

    expression = _expression_frame(expression, atlas_info)
    _write_expression_archive(expression, expression_out)
    counts.to_csv(counts_out, index=True)
    if Path(atlas_info).exists():
        pd.read_csv(atlas_info).to_csv(atlas_info_out, index=False)
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
            },
            indent=2,
        )
    )
    return {
        "expression": expression_out,
        "counts": counts_out,
        "report": report_out,
        "provenance": provenance_out,
        "labels": atlas_info_out,
    }
