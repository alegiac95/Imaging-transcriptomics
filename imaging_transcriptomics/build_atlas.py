from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

import pandas as pd

from .atlas_registry import get_atlas


MirrorMode = Literal[None, "bidirectional", "leftright", "rightleft"]


def _import_abagen():
    try:
        import abagen
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "abagen is required to build atlas expression assets. Install imaging-transcriptomics[maps]."
        ) from exc
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
) -> dict[str, Path]:
    abagen = _import_abagen()
    spec = get_atlas(atlas)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    atlas_image = atlas_image or spec.volume_1mm_path
    atlas_info = atlas_info or spec.labels_path
    if atlas_image is None or atlas_info is None:
        raise ValueError(
            f"Atlas '{spec.id}' does not have enough packaged metadata to build expression assets automatically. "
            "Provide atlas_image and atlas_info explicitly."
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

    expression_out = output_path / f"atlas-{spec.id}_gene_expression_data.csv"
    counts_out = output_path / f"atlas-{spec.id}_sample_counts.csv"
    report_out = output_path / "README.txt"
    provenance_out = output_path / "provenance.json"
    atlas_info_out = output_path / f"atlas-{spec.id}_labels.csv"

    expression = expression.reset_index().rename(columns={expression.index.name or "index": "id"})
    if "label" not in expression.columns:
        labels = pd.read_csv(atlas_info)
        if "Unnamed: 0" in labels.columns:
            labels = labels.drop(columns=["Unnamed: 0"])
        expression = labels[["id", "label"]].merge(expression, on="id", how="right")
    expression.to_csv(expression_out, index=False)
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
