from __future__ import annotations

import argparse
import shutil
import ast
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from netneurotools import datasets as nnt_datasets
from neuromaps.datasets import fetch_atlas as fetch_nm_atlas
from nilearn.datasets import fetch_atlas_destrieux_2009, fetch_atlas_schaefer_2018
from nilearn.image import resample_to_img

import imaging_transcriptomics as imt


def _copy_file(src, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    print(f"copied {src} -> {dst}", flush=True)


def _write_labels(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"wrote {path} ({len(df)} rows)", flush=True)


def _prepare_schaefer(
    atlas_root: Path,
    tmp_nilearn: Path,
    tmp_nnt: Path,
) -> None:
    surfaces = nnt_datasets.fetch_schaefer2018(version="fsaverage", data_dir=tmp_nnt, verbose=1)
    for n_rois in (200, 400):
        outdir = atlas_root / f"Schaefer_{n_rois}"
        outdir.mkdir(parents=True, exist_ok=True)
        sch_1mm = fetch_atlas_schaefer_2018(
            n_rois=n_rois,
            resolution_mm=1,
            data_dir=tmp_nilearn,
            verbose=1,
        )
        sch_2mm = fetch_atlas_schaefer_2018(
            n_rois=n_rois,
            resolution_mm=2,
            data_dir=tmp_nilearn,
            verbose=1,
        )
        surface = surfaces[f"{n_rois}Parcels7Networks"]
        _copy_file(sch_1mm["maps"], outdir / f"atlas-Schaefer_{n_rois}_1mm.nii.gz")
        _copy_file(sch_2mm["maps"], outdir / f"atlas-Schaefer_{n_rois}_2mm.nii.gz")
        _copy_file(surface.L, outdir / f"atlas-Schaefer_{n_rois}_lh_aparc.annot")
        _copy_file(surface.R, outdir / f"atlas-Schaefer_{n_rois}_rh_aparc.annot")

        raw_labels = list(sch_1mm["labels"])
        if len(raw_labels) == n_rois + 1:
            raw_labels = raw_labels[1:]
        labels = []
        for idx, label in enumerate(raw_labels, start=1):
            name = label.decode() if isinstance(label, bytes) else str(label)
            hemisphere = "L" if "_LH_" in name else "R" if "_RH_" in name else "B"
            labels.append(
                {
                    "id": idx,
                    "label": name,
                    "hemisphere": hemisphere,
                    "structure": "cortex",
                }
            )
        _write_labels(pd.DataFrame(labels), outdir / f"atlas-schaefer-{n_rois}_labels.csv")


def _prepare_destrieux(atlas_root: Path, tmp_nilearn: Path) -> None:
    def _parse_label(raw) -> tuple[str, str]:
        value = raw
        if isinstance(raw, str):
            try:
                parsed = ast.literal_eval(raw)
            except (SyntaxError, ValueError):
                parsed = raw
            value = parsed
        if isinstance(value, tuple) and len(value) >= 2:
            text = str(value[1])
        else:
            text = str(value)
        hemisphere = "L" if text.startswith("L ") else "R" if text.startswith("R ") else "B"
        label = text[2:] if hemisphere in {"L", "R"} and len(text) > 2 else text
        return label, hemisphere

    outdir = atlas_root / "Destrieux"
    outdir.mkdir(parents=True, exist_ok=True)
    destrieux = fetch_atlas_destrieux_2009(lateralized=True, data_dir=tmp_nilearn, verbose=1)
    _copy_file(destrieux["maps"], outdir / "atlas-Destrieux_1mm.nii.gz")

    ref_2mm = nib.load(str(atlas_root / "DK" / "atlas-DK_2mm.nii.gz"))
    destrieux_img = nib.load(str(destrieux["maps"]))
    destrieux_2mm = resample_to_img(
        destrieux_img,
        ref_2mm,
        interpolation="nearest",
        force_resample=True,
    )
    destrieux_2mm_path = outdir / "atlas-Destrieux_2mm.nii.gz"
    nib.save(destrieux_2mm, str(destrieux_2mm_path))
    print(f"wrote {destrieux_2mm_path}", flush=True)

    present_ids = set(int(x) for x in np.unique(destrieux_img.get_fdata()) if x > 0)
    if "lut" in destrieux:
        lut = destrieux["lut"].copy()
        label_col = "name" if "name" in lut.columns else lut.columns[-1]
        id_col = "index" if "index" in lut.columns else lut.columns[0]
    else:
        labels_raw = list(destrieux["labels"])
        if labels_raw and str(labels_raw[0]).lower() == "background":
            labels_raw = labels_raw[1:]
            ids = range(1, len(labels_raw) + 1)
        else:
            ids = range(len(labels_raw))
        lut = pd.DataFrame({"id": list(ids), "label": [str(label) for label in labels_raw]})
        label_col = "label"
        id_col = "id"
    lut = lut.loc[lut[id_col].isin(sorted(present_ids))].copy()
    parsed = lut[label_col].apply(_parse_label)
    labels = pd.DataFrame(
        {
            "id": lut[id_col].astype(int).to_numpy(),
            "label": [label for label, _ in parsed],
            "hemisphere": [hemisphere for _, hemisphere in parsed],
        }
    )
    labels["structure"] = "cortex"
    _write_labels(labels.sort_values("id").reset_index(drop=True), outdir / "atlas-destrieux_labels.csv")


def _prepare_glasser(atlas_root: Path, tmp_nnt: Path, tmp_nm: Path) -> None:
    outdir = atlas_root / "Glasser_360"
    outdir.mkdir(parents=True, exist_ok=True)

    mmp = nnt_datasets.fetch_mmpall(data_dir=tmp_nnt, verbose=1)
    fslr = fetch_nm_atlas("fslr", "32k", data_dir=tmp_nm, verbose=1)

    _copy_file(mmp.L, outdir / "atlas-Glasser_360_lh.label.gii")
    _copy_file(mmp.R, outdir / "atlas-Glasser_360_rh.label.gii")
    _copy_file(fslr["midthickness"].L, outdir / "atlas-Glasser_360_lh.surf.gii")
    _copy_file(fslr["midthickness"].R, outdir / "atlas-Glasser_360_rh.surf.gii")

    label_img = nib.load(str(mmp.L))
    label_table = label_img.labeltable.get_labels_as_dict()
    labels = []
    for idx in sorted(k for k in label_table if int(k) > 0):
        name = str(label_table[idx])
        hemisphere = "L" if name.startswith("L_") else "R" if name.startswith("R_") else "B"
        labels.append(
            {
                "id": int(idx),
                "label": name,
                "hemisphere": hemisphere,
                "structure": "cortex",
            }
        )
    _write_labels(pd.DataFrame(labels), outdir / "atlas-glasser-360_labels.csv")


def _build_expression(
    atlas_root: Path,
    abagen_cache: Path,
    n_proc: int,
    atlas_ids: set[str] | None = None,
) -> None:
    jobs = [
        ("dk", atlas_root / "DK", None, None, None, None),
        ("schaefer-100", atlas_root / "Schaefer_100", None, None, None, None),
        ("schaefer-200", atlas_root / "Schaefer_200", None, None, None, None),
        ("schaefer-400", atlas_root / "Schaefer_400", None, None, None, None),
        ("destrieux", atlas_root / "Destrieux", None, None, None, None),
        (
            "glasser-360",
            atlas_root / "Glasser_360",
            (
                atlas_root / "Glasser_360" / "atlas-Glasser_360_lh.label.gii",
                atlas_root / "Glasser_360" / "atlas-Glasser_360_rh.label.gii",
            ),
            atlas_root / "Glasser_360" / "atlas-glasser-360_labels.csv",
            (
                atlas_root / "Glasser_360" / "atlas-Glasser_360_lh.surf.gii",
                atlas_root / "Glasser_360" / "atlas-Glasser_360_rh.surf.gii",
            ),
            "fslr",
        ),
    ]

    for atlas, outdir, atlas_image, atlas_info, geometry, space in jobs:
        if atlas_ids is not None and atlas not in atlas_ids:
            continue
        print(f"\n=== BUILD {atlas} ===", flush=True)
        kwargs = dict(
            atlas=atlas,
            output_dir=outdir,
            lr_mirror="leftright",
            donors="all",
            n_proc=n_proc,
            data_dir=abagen_cache,
        )
        if atlas_image is not None:
            kwargs["atlas_image"] = atlas_image
        if atlas_info is not None:
            kwargs["atlas_info"] = atlas_info
        if geometry is not None:
            kwargs["geometry"] = geometry
        if space is not None:
            kwargs["space"] = space

        outputs = imt.build_expression_assets(**kwargs)
        for key, value in outputs.items():
            print(f"{atlas} {key}: {value}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch and build all preset atlas assets.")
    parser.add_argument(
        "--n-proc",
        type=int,
        default=4,
        help="Number of worker processes to use inside abagen.",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Repository root containing imaging_transcriptomics/.",
    )
    parser.add_argument(
        "--atlases",
        nargs="*",
        default=None,
        help="Optional subset of atlas ids to prepare/build, e.g. dk schaefer-200 glasser-360.",
    )
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    atlas_root = repo_root / "imaging_transcriptomics" / "data" / "atlases"
    tmp_nilearn = repo_root / "tmp_atlas_fetch"
    tmp_nnt = repo_root / "tmp_atlas_fetch_nnt"
    tmp_nm = repo_root / "tmp_atlas_fetch_neuromaps"
    abagen_cache = repo_root / ".cache" / "abagen-data"

    for path in (atlas_root, tmp_nilearn, tmp_nnt, tmp_nm, abagen_cache):
        path.mkdir(parents=True, exist_ok=True)

    requested = set(args.atlases or ["dk", "schaefer-100", "schaefer-200", "schaefer-400", "destrieux", "glasser-360"])
    if requested.intersection({"schaefer-200", "schaefer-400"}):
        _prepare_schaefer(atlas_root, tmp_nilearn, tmp_nnt)
    if "destrieux" in requested:
        _prepare_destrieux(atlas_root, tmp_nilearn)
    if "glasser-360" in requested:
        _prepare_glasser(atlas_root, tmp_nnt, tmp_nm)
    if requested.intersection({"dk", "schaefer-100", "schaefer-200", "schaefer-400", "destrieux", "glasser-360"}):
        _build_expression(atlas_root, abagen_cache, n_proc=args.n_proc, atlas_ids=requested)
    print("\nall preset atlas builds complete", flush=True)


if __name__ == "__main__":
    main()
