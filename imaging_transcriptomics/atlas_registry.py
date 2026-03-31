from __future__ import annotations

from pathlib import Path

import pandas as pd

from .models import AtlasSpec


ATLAS_DATA_DIR = Path(__file__).resolve().parent / "data" / "atlases"
ATLAS_SHARED_DIR = ATLAS_DATA_DIR / "shared"


def _atlas_path(*parts: str) -> Path:
    return ATLAS_DATA_DIR.joinpath(*parts)


def _shared_path(*parts: str) -> Path:
    return ATLAS_SHARED_DIR.joinpath(*parts)


ATLAS_REGISTRY: dict[str, AtlasSpec] = {
    "dk": AtlasSpec(
        id="dk",
        label="Desikan-Killiany",
        description="Legacy-compatible Desikan-Killiany atlas with cortical parcels and the packaged aseg-derived subcortical add-on.",
        family="Desikan-Killiany",
        volumetric_space="MNI152",
        supported_spaces=("MNI152", "fsaverage"),
        packaged=True,
        buildable=True,
        default_hemisphere="left",
        has_subcortex=True,
        n_regions_left=41,
        n_regions_both=83,
        labels_path=_atlas_path("DK", "atlas-DK_labels.csv"),
        expression_path=_atlas_path("DK", "atlas-DK_gene_expression_data.npz"),
        gene_labels_path=_shared_path("genes-ahba-15677.npy"),
        volume_1mm_path=_atlas_path("DK", "atlas-DK_1mm.nii.gz"),
        volume_2mm_path=_atlas_path("DK", "atlas-DK_2mm.nii.gz"),
        lh_annot_path=_atlas_path("DK", "atlas-DK_fsa5_lh_aparc.annot"),
        rh_annot_path=_atlas_path("DK", "atlas-DK_fsa5_rh_aparc.annot"),
        surface_space="fsaverage",
        surface_density="10k",
        notes="The packaged atlas includes the aseg subcortical add-on. Use regions='default' or 'cort' for cortex only, and regions='all' or 'cort+sub' for whole-brain coverage.",
    ),
    "schaefer-100": AtlasSpec(
        id="schaefer-100",
        label="Schaefer 100",
        description="7-network Schaefer atlas with packaged cortical parcels plus the same aseg-derived subcortical add-on used in dk.",
        family="Schaefer",
        volumetric_space="MNI152",
        supported_spaces=("MNI152", "fsaverage"),
        packaged=True,
        buildable=True,
        default_hemisphere="left",
        has_subcortex=True,
        n_regions_left=57,
        n_regions_both=115,
        labels_path=_atlas_path("Schaefer_100_Aseg", "atlas-schaefer-100-aseg_labels.csv"),
        expression_path=_atlas_path(
            "Schaefer_100_Aseg", "atlas-schaefer-100-aseg_gene_expression_data.npz"
        ),
        gene_labels_path=_shared_path("genes-ahba-15677.npy"),
        volume_1mm_path=_atlas_path("Schaefer_100_Aseg", "atlas-Schaefer_100_1mm.nii.gz"),
        volume_2mm_path=_atlas_path("Schaefer_100_Aseg", "atlas-Schaefer_100_2mm.nii.gz"),
        lh_annot_path=_atlas_path("Schaefer_100_Aseg", "atlas-Schaefer_100_lh_aparc.annot"),
        rh_annot_path=_atlas_path("Schaefer_100_Aseg", "atlas-Schaefer_100_rh_aparc.annot"),
        surface_space="fsaverage",
        surface_density="10k",
        notes="Use regions='default' or 'cort' for the original 100 cortical parcels, and regions='all' or 'cort+sub' to append the shared aseg subcortical add-on.",
    ),
    "schaefer-200": AtlasSpec(
        id="schaefer-200",
        label="Schaefer 200",
        description="Higher-resolution Schaefer atlas with packaged cortical parcels plus the same aseg-derived subcortical add-on used in dk.",
        family="Schaefer",
        volumetric_space="MNI152",
        supported_spaces=("MNI152", "fsaverage", "fsLR"),
        packaged=True,
        buildable=True,
        default_hemisphere="left",
        has_subcortex=True,
        n_regions_left=107,
        n_regions_both=215,
        labels_path=_atlas_path("Schaefer_200_Aseg", "atlas-schaefer-200-aseg_labels.csv"),
        expression_path=_atlas_path(
            "Schaefer_200_Aseg", "atlas-schaefer-200-aseg_gene_expression_data.npz"
        ),
        gene_labels_path=_shared_path("genes-ahba-15677.npy"),
        volume_1mm_path=_atlas_path("Schaefer_200_Aseg", "atlas-Schaefer_200_1mm.nii.gz"),
        volume_2mm_path=_atlas_path("Schaefer_200_Aseg", "atlas-Schaefer_200_2mm.nii.gz"),
        lh_annot_path=_atlas_path("Schaefer_200_Aseg", "atlas-Schaefer_200_lh_aparc.annot"),
        rh_annot_path=_atlas_path("Schaefer_200_Aseg", "atlas-Schaefer_200_rh_aparc.annot"),
        surface_space="fsaverage",
        surface_density="10k",
        notes="Use regions='default' or 'cort' for the original 200 cortical parcels, and regions='all' or 'cort+sub' to append the shared aseg subcortical add-on.",
    ),
    "schaefer-400": AtlasSpec(
        id="schaefer-400",
        label="Schaefer 400",
        description="Higher-resolution Schaefer atlas with packaged cortical parcels plus the same aseg-derived subcortical add-on used in dk.",
        family="Schaefer",
        volumetric_space="MNI152",
        supported_spaces=("MNI152", "fsaverage", "fsLR"),
        packaged=True,
        buildable=True,
        default_hemisphere="left",
        has_subcortex=True,
        n_regions_left=207,
        n_regions_both=415,
        labels_path=_atlas_path("Schaefer_400_Aseg", "atlas-schaefer-400-aseg_labels.csv"),
        expression_path=_atlas_path(
            "Schaefer_400_Aseg", "atlas-schaefer-400-aseg_gene_expression_data.npz"
        ),
        gene_labels_path=_shared_path("genes-ahba-15675.npy"),
        volume_1mm_path=_atlas_path("Schaefer_400_Aseg", "atlas-Schaefer_400_1mm.nii.gz"),
        volume_2mm_path=_atlas_path("Schaefer_400_Aseg", "atlas-Schaefer_400_2mm.nii.gz"),
        lh_annot_path=_atlas_path("Schaefer_400_Aseg", "atlas-Schaefer_400_lh_aparc.annot"),
        rh_annot_path=_atlas_path("Schaefer_400_Aseg", "atlas-Schaefer_400_rh_aparc.annot"),
        surface_space="fsaverage",
        surface_density="10k",
        notes="Use regions='default' or 'cort' for the original 400 cortical parcels, and regions='all' or 'cort+sub' to append the shared aseg subcortical add-on.",
    ),
    "destrieux": AtlasSpec(
        id="destrieux",
        label="Destrieux",
        description="Destrieux atlas with packaged cortical parcels plus the same aseg-derived subcortical add-on used in dk.",
        family="Destrieux",
        volumetric_space="MNI152",
        supported_spaces=("MNI152", "fsaverage"),
        packaged=True,
        buildable=True,
        default_hemisphere="left",
        has_subcortex=True,
        n_regions_left=81,
        n_regions_both=163,
        labels_path=_atlas_path("Destrieux_Aseg", "atlas-destrieux-aseg_labels.csv"),
        expression_path=_atlas_path("Destrieux_Aseg", "atlas-destrieux-aseg_gene_expression_data.npz"),
        gene_labels_path=_shared_path("genes-ahba-15675.npy"),
        volume_1mm_path=_atlas_path("Destrieux_Aseg", "atlas-Destrieux_1mm.nii.gz"),
        volume_2mm_path=_atlas_path("Destrieux_Aseg", "atlas-Destrieux_2mm.nii.gz"),
        surface_space="fsaverage",
        surface_density="10k",
        notes="Use regions='default' or 'cort' for the original Destrieux cortex, and regions='all' or 'cort+sub' to append the shared aseg subcortical add-on.",
    ),
    "glasser-360": AtlasSpec(
        id="glasser-360",
        label="Glasser 360",
        description="Surface-oriented Glasser atlas with packaged cortical parcels plus the same aseg-derived subcortical add-on used in dk.",
        family="Glasser",
        volumetric_space="MNI152",
        supported_spaces=("fsLR", "fsaverage"),
        packaged=True,
        buildable=True,
        default_hemisphere="left",
        has_subcortex=True,
        n_regions_left=187,
        n_regions_both=375,
        labels_path=_atlas_path("Glasser_360_Aseg", "atlas-glasser-360-aseg_labels.csv"),
        expression_path=_atlas_path("Glasser_360_Aseg", "atlas-glasser-360-aseg_gene_expression_data.npz"),
        gene_labels_path=_shared_path("genes-ahba-15675.npy"),
        lh_surface_path=_atlas_path("Glasser_360_Aseg", "atlas-Glasser_360_lh.label.gii"),
        rh_surface_path=_atlas_path("Glasser_360_Aseg", "atlas-Glasser_360_rh.label.gii"),
        geometry_lh_path=_atlas_path("Glasser_360_Aseg", "atlas-Glasser_360_lh.surf.gii"),
        geometry_rh_path=_atlas_path("Glasser_360_Aseg", "atlas-Glasser_360_rh.surf.gii"),
        surface_space="fslr",
        surface_density="32k",
        notes="Use regions='default' or 'cort' for the original 360 cortical parcels. Regions='all' or 'cort+sub' append the shared aseg subcortical add-on; this packaged hybrid remains surface-first.",
    ),
}


LEGACY_ALIASES = {
    "DK": "dk",
    "dk": "dk",
    "Schaefer_100": "schaefer-100",
    "schaefer_100": "schaefer-100",
    "schaefer-100": "schaefer-100",
    "Schaefer_200": "schaefer-200",
    "schaefer_200": "schaefer-200",
    "schaefer-200": "schaefer-200",
    "Schaefer_400": "schaefer-400",
    "schaefer_400": "schaefer-400",
    "schaefer-400": "schaefer-400",
    "destrieux": "destrieux",
    "glasser360": "glasser-360",
    "glasser-360": "glasser-360",
    "Schaefer_100_Aseg": "schaefer-100",
    "schaefer_100_aseg": "schaefer-100",
    "schaefer-100-aseg": "schaefer-100",
    "Schaefer_200_Aseg": "schaefer-200",
    "schaefer_200_aseg": "schaefer-200",
    "schaefer-200-aseg": "schaefer-200",
    "Schaefer_400_Aseg": "schaefer-400",
    "schaefer_400_aseg": "schaefer-400",
    "schaefer-400-aseg": "schaefer-400",
    "destrieux_aseg": "destrieux",
    "destrieux-aseg": "destrieux",
    "glasser_360_aseg": "glasser-360",
    "glasser360aseg": "glasser-360",
    "glasser-360-aseg": "glasser-360",
}


def normalize_atlas_id(atlas: str) -> str:
    try:
        return LEGACY_ALIASES[atlas]
    except KeyError as exc:
        valid = ", ".join(sorted(ATLAS_REGISTRY))
        raise ValueError(f"Unknown atlas '{atlas}'. Available atlases: {valid}.") from exc


def get_atlas(atlas: str) -> AtlasSpec:
    return ATLAS_REGISTRY[normalize_atlas_id(atlas)]


def list_atlases(packaged_only: bool = False) -> list[AtlasSpec]:
    atlases = list(ATLAS_REGISTRY.values())
    if packaged_only:
        atlases = [atlas for atlas in atlases if atlas.packaged]
    return atlases


def describe_atlas(atlas: str) -> dict[str, object]:
    spec = get_atlas(atlas)
    return {
        "id": spec.id,
        "label": spec.label,
        "description": spec.description,
        "family": spec.family,
        "packaged": spec.packaged,
        "buildable": spec.buildable,
        "default_hemisphere": spec.default_hemisphere,
        "has_subcortex": spec.has_subcortex,
        "n_regions_left": spec.n_regions_left,
        "n_regions_both": spec.n_regions_both,
        "supported_spaces": spec.supported_spaces,
        "notes": spec.notes,
    }


def atlas_table(packaged_only: bool = False) -> pd.DataFrame:
    return pd.DataFrame.from_records(
        describe_atlas(atlas.id) for atlas in list_atlases(packaged_only=packaged_only)
    )
