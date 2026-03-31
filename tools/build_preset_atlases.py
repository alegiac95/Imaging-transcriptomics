from __future__ import annotations

import argparse
import ast
import json
import shutil
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from netneurotools import datasets as nnt_datasets
from neuromaps.datasets import fetch_atlas as fetch_nm_atlas
from nilearn.datasets import fetch_atlas_destrieux_2009, fetch_atlas_schaefer_2018
from nilearn.image import resample_to_img

import imaging_transcriptomics as imt


HYBRID_ATLAS_CONFIG = {
    "schaefer-100-aseg": {"base_id": "schaefer-100", "outdir_name": "Schaefer_100_Aseg"},
    "schaefer-200-aseg": {"base_id": "schaefer-200", "outdir_name": "Schaefer_200_Aseg"},
    "schaefer-400-aseg": {"base_id": "schaefer-400", "outdir_name": "Schaefer_400_Aseg"},
    "destrieux-aseg": {"base_id": "destrieux", "outdir_name": "Destrieux_Aseg"},
    "glasser-360-aseg": {"base_id": "glasser-360", "outdir_name": "Glasser_360_Aseg"},
}


def _copy_file(src, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    print(f"copied {src} -> {dst}", flush=True)


def _write_labels(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"wrote {path} ({len(df)} rows)", flush=True)


def _read_labels(path: Path) -> pd.DataFrame:
    """Load one atlas labels table into a normalized DataFrame."""

    labels = pd.read_csv(path)
    if "Unnamed: 0" in labels.columns:
        labels = labels.drop(columns=["Unnamed: 0"])
    return labels


def _read_expression_archive(path: Path) -> np.ndarray:
    """Load one packaged NPZ expression matrix as float32 values."""

    with np.load(path, allow_pickle=False) as archive:
        return archive["values"].astype(np.float32, copy=False)


def _write_expression_archive(values: np.ndarray, path: Path) -> None:
    """Write a packaged NPZ expression archive."""

    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, values=np.asarray(values, dtype=np.float32, copy=False))
    print(f"wrote {path} ({values.shape[0]} x {values.shape[1]})", flush=True)


def _copy_optional(src: Path | None, dst_dir: Path) -> Path | None:
    """Copy one asset into ``dst_dir`` when the source exists."""

    if src is None:
        return None
    dst = dst_dir / src.name
    _copy_file(src, dst)
    return dst


def _merge_aseg_volume(
    *,
    base_volume_path: Path,
    dk_volume_path: Path,
    id_mapping: dict[int, int],
    output_path: Path,
) -> None:
    """Append dk/aseg subcortex labels to a cortical atlas volume."""

    base_img = nib.load(str(base_volume_path))
    base_data = np.asarray(base_img.dataobj, dtype=np.int32).copy()
    dk_img = nib.load(str(dk_volume_path))
    if dk_img.shape != base_img.shape or not np.allclose(dk_img.affine, base_img.affine):
        dk_img = resample_to_img(
            dk_img,
            base_img,
            interpolation="nearest",
            force_resample=True,
        )
    dk_data = np.asarray(dk_img.dataobj, dtype=np.int32)
    for source_id, target_id in id_mapping.items():
        base_data[dk_data == source_id] = target_id
    merged = nib.Nifti1Image(base_data, base_img.affine, header=base_img.header)
    merged.header.set_data_dtype(np.int32)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(merged, str(output_path))
    print(f"wrote {output_path}", flush=True)


def _match_sample_counts(base_path: Path | None, dk_path: Path | None) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """Load optional sample-count tables when they are available."""

    if base_path is None or dk_path is None or not base_path.exists() or not dk_path.exists():
        return None, None
    base = pd.read_csv(base_path, index_col=0).reset_index().rename(columns={"index": "id"})
    dk = pd.read_csv(dk_path, index_col=0).reset_index().rename(columns={"index": "id"})
    return base, dk


def _sample_counts_path(spec) -> Path | None:
    """Return the packaged sample-count table next to one atlas, when present."""

    if spec.expression_path is None:
        return None
    matches = sorted(spec.expression_path.parent.glob("*sample_counts.csv"))
    return matches[0] if matches else None


def _write_hybrid_metadata(
    *,
    output_dir: Path,
    hybrid_id: str,
    base_id: str,
    dk_sub_labels: pd.DataFrame,
    volume_packaged: bool,
    gene_count: int,
    base_readme: str | None,
) -> None:
    """Write a short README and provenance file for one synthetic hybrid atlas."""

    lines = [
        f"Hybrid preset: {hybrid_id}",
        "",
        f"This packaged atlas uses the cortical parcels from '{base_id}' and appends the same aseg-derived subcortical add-on used in 'dk'.",
        "The cortical expression matrix is copied from the packaged base atlas.",
        "The subcortical expression rows are copied from the packaged dk atlas after aligning them to the base atlas gene ordering.",
        f"Packaged genes: {gene_count}",
        f"Transferred aseg parcels: {len(dk_sub_labels)}",
        f"Volumetric merged atlas available: {'yes' if volume_packaged else 'no'}",
    ]
    if base_readme:
        lines.extend(
            [
                "",
                "Base atlas notes:",
                base_readme.strip(),
            ]
        )
    (output_dir / "README.txt").write_text("\n".join(lines) + "\n")
    (output_dir / "provenance.json").write_text(
        json.dumps(
            {
                "atlas": hybrid_id,
                "base_atlas": base_id,
                "subcortex_source": "dk",
                "subcortical_labels": dk_sub_labels["label"].astype(str).tolist(),
                "subcortical_hemispheres": dk_sub_labels["hemisphere"].astype(str).tolist(),
                "volume_packaged": volume_packaged,
                "expression_strategy": "base cortical matrix plus dk aseg rows aligned to base gene labels",
            },
            indent=2,
        )
    )


def _build_aseg_hybrid(atlas_root: Path, hybrid_id: str) -> None:
    """Synthesize a packaged ``*-aseg`` atlas from local packaged assets."""

    config = HYBRID_ATLAS_CONFIG[hybrid_id]
    base_spec = imt.get_atlas(config["base_id"])
    hybrid_spec = imt.get_atlas(hybrid_id)
    dk_spec = imt.get_atlas("dk")

    if base_spec.labels_path is None or base_spec.expression_path is None or base_spec.gene_labels_path is None:
        raise RuntimeError(f"Base atlas '{base_spec.id}' is missing packaged assets required for hybrid synthesis.")
    if dk_spec.labels_path is None or dk_spec.expression_path is None or dk_spec.gene_labels_path is None:
        raise RuntimeError("The packaged dk atlas is missing the aseg assets required for hybrid synthesis.")

    output_dir = atlas_root / config["outdir_name"]
    output_dir.mkdir(parents=True, exist_ok=True)

    base_labels = _read_labels(base_spec.labels_path)
    dk_labels = _read_labels(dk_spec.labels_path)
    dk_sub_labels = dk_labels.loc[dk_labels["structure"].astype(str) != "cortex"].copy().reset_index(drop=True)
    next_ids = np.arange(int(base_labels["id"].max()) + 1, int(base_labels["id"].max()) + 1 + len(dk_sub_labels))
    hybrid_sub_labels = dk_sub_labels.copy()
    hybrid_sub_labels["id"] = next_ids
    hybrid_labels = pd.concat([base_labels, hybrid_sub_labels], ignore_index=True)
    _write_labels(hybrid_labels, hybrid_spec.labels_path)

    base_genes = np.load(base_spec.gene_labels_path, allow_pickle=False).astype(str, copy=False)
    dk_genes = np.load(dk_spec.gene_labels_path, allow_pickle=False).astype(str, copy=False)
    base_values = _read_expression_archive(base_spec.expression_path)
    dk_values = _read_expression_archive(dk_spec.expression_path)
    dk_gene_index = {gene: idx for idx, gene in enumerate(dk_genes)}
    try:
        dk_column_index = np.asarray([dk_gene_index[gene] for gene in base_genes], dtype=np.int32)
    except KeyError as exc:
        raise RuntimeError(
            f"Unable to align dk aseg expression rows to '{base_spec.id}' because gene '{exc.args[0]}' is missing from the dk shared gene list."
        ) from exc
    dk_sub_row_index = dk_labels.index[dk_labels["structure"].astype(str) != "cortex"].to_numpy(dtype=np.int32, copy=False)
    hybrid_sub_values = dk_values[np.ix_(dk_sub_row_index, dk_column_index)]
    hybrid_values = np.vstack([base_values, hybrid_sub_values]).astype(np.float32, copy=False)
    _write_expression_archive(hybrid_values, hybrid_spec.expression_path)

    base_counts, dk_counts = _match_sample_counts(_sample_counts_path(base_spec), _sample_counts_path(dk_spec))
    if base_counts is not None and dk_counts is not None:
        source_to_target = dict(zip(dk_sub_labels["id"].astype(int), next_ids, strict=True))
        hybrid_sub_counts = dk_counts.loc[dk_counts["id"].astype(int).isin(source_to_target)].copy()
        hybrid_sub_counts["id"] = hybrid_sub_counts["id"].astype(int).map(source_to_target)
        pd.concat([base_counts, hybrid_sub_counts], ignore_index=True).to_csv(
            output_dir / f"atlas-{hybrid_id}_sample_counts.csv",
            index=False,
        )

    volume_packaged = False
    if base_spec.volume_1mm_path is not None and base_spec.volume_2mm_path is not None:
        source_to_target = dict(zip(dk_sub_labels["id"].astype(int), next_ids, strict=True))
        _merge_aseg_volume(
            base_volume_path=base_spec.volume_1mm_path,
            dk_volume_path=dk_spec.volume_1mm_path,
            id_mapping=source_to_target,
            output_path=hybrid_spec.volume_1mm_path,
        )
        _merge_aseg_volume(
            base_volume_path=base_spec.volume_2mm_path,
            dk_volume_path=dk_spec.volume_2mm_path,
            id_mapping=source_to_target,
            output_path=hybrid_spec.volume_2mm_path,
        )
        volume_packaged = True

    for source in (
        base_spec.lh_annot_path,
        base_spec.rh_annot_path,
        base_spec.lh_surface_path,
        base_spec.rh_surface_path,
        base_spec.geometry_lh_path,
        base_spec.geometry_rh_path,
    ):
        _copy_optional(source, output_dir)

    base_readme_path = base_spec.expression_path.parent / "README.txt"
    base_readme = base_readme_path.read_text() if base_readme_path.exists() else None
    _write_hybrid_metadata(
        output_dir=output_dir,
        hybrid_id=hybrid_id,
        base_id=base_spec.id,
        dk_sub_labels=dk_sub_labels,
        volume_packaged=volume_packaged,
        gene_count=int(base_genes.shape[0]),
        base_readme=base_readme,
    )


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


def _build_hybrid_atlases(atlas_root: Path, atlas_ids: set[str] | None = None) -> None:
    """Generate all requested ``*-aseg`` atlas assets from packaged sources."""

    for hybrid_id in HYBRID_ATLAS_CONFIG:
        if atlas_ids is not None and hybrid_id not in atlas_ids:
            continue
        print(f"\n=== SYNTHESIZE {hybrid_id} ===", flush=True)
        _build_aseg_hybrid(atlas_root, hybrid_id)


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

    requested = set(
        args.atlases
        or [
            "dk",
            "schaefer-100",
            "schaefer-200",
            "schaefer-400",
            "destrieux",
            "glasser-360",
            *HYBRID_ATLAS_CONFIG,
        ]
    )
    if requested.intersection({"schaefer-200", "schaefer-400"}):
        _prepare_schaefer(atlas_root, tmp_nilearn, tmp_nnt)
    if "destrieux" in requested:
        _prepare_destrieux(atlas_root, tmp_nilearn)
    if "glasser-360" in requested:
        _prepare_glasser(atlas_root, tmp_nnt, tmp_nm)
    if requested.intersection({"dk", "schaefer-100", "schaefer-200", "schaefer-400", "destrieux", "glasser-360"}):
        _build_expression(atlas_root, abagen_cache, n_proc=args.n_proc, atlas_ids=requested)
    if requested.intersection(HYBRID_ATLAS_CONFIG):
        _build_hybrid_atlases(atlas_root, atlas_ids=requested)
    print("\nall preset atlas builds complete", flush=True)


if __name__ == "__main__":
    main()
