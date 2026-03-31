from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from imaging_transcriptomics import extract_scan_data
from imaging_transcriptomics.exceptions import InputDataError
from imaging_transcriptomics.inputs.common import (
    load_tabular_vector,
    make_extracted_scan,
    recognized_suffix,
    surface_inputs,
)
import imaging_transcriptomics.inputs.neuromaps as neuromaps_inputs
import imaging_transcriptomics.inputs.volume as volume_inputs
import imaging_transcriptomics.scan as scan_module


@pytest.mark.parametrize(
    ("filename", "expected"),
    [
        ("subject.native.T1w.nii.gz", ".nii.gz"),
        ("left.func.gii", ".func.gii"),
        ("right.shape.gii", ".shape.gii"),
        ("regional_vector.tsv", ".tsv"),
        ("unknown.ext", ".ext"),
    ],
)
def test_recognized_suffix_handles_multidot_and_known_variants(filename: str, expected: str):
    assert recognized_suffix(Path(filename)) == expected


def test_load_tabular_vector_accepts_single_column_and_rejects_matrices(tmp_path: Path):
    vector_path = tmp_path / "vector.txt"
    matrix_path = tmp_path / "matrix.csv"
    np.savetxt(vector_path, np.array([1.0, 2.0, 3.0], dtype=float))
    np.savetxt(matrix_path, np.array([[1.0, 2.0], [3.0, 4.0]], dtype=float), delimiter=",")

    loaded = load_tabular_vector(vector_path)

    np.testing.assert_allclose(loaded, np.array([1.0, 2.0, 3.0], dtype=float))
    with pytest.raises(InputDataError, match="one-column vector"):
        load_tabular_vector(matrix_path)


def test_surface_inputs_normalizes_supported_call_styles():
    assert surface_inputs(("lh.func.gii", "rh.func.gii"), None) == ("lh.func.gii", "rh.func.gii")
    assert surface_inputs("lh.func.gii", "rh.func.gii") == ("lh.func.gii", "rh.func.gii")
    assert surface_inputs("lh.func.gii", None) is None


def test_make_extracted_scan_rejects_wrong_vector_length(dk_left_selection):
    with pytest.raises(InputDataError, match="Expected 41 regional values"):
        make_extracted_scan(
            dk_left_selection,
            np.arange(40, dtype=float),
            source="array",
            source_space=None,
            source_kind="vector",
        )


def test_extract_scan_data_dispatches_text_vectors(tmp_path: Path, dk_left_selection):
    values = np.linspace(-1.0, 1.0, dk_left_selection.n_regions, dtype=float)
    path = tmp_path / "regional_values.txt"
    np.savetxt(path, values)

    extracted = extract_scan_data(path, atlas="dk", hemisphere="left", regions="all")

    np.testing.assert_allclose(extracted.values, values)
    assert extracted.source_kind == "vector"


def test_extract_scan_data_dispatches_surface_pairs(monkeypatch, dk_left_selection):
    expected = np.linspace(-1.0, 1.0, dk_left_selection.n_regions, dtype=float)
    seen = {}

    def fake_parcellate(data, selection, *, source_space, hemi=None):
        seen["data"] = data
        seen["selection"] = selection
        seen["source_space"] = source_space
        seen["hemi"] = hemi
        return expected

    monkeypatch.setattr(scan_module, "parcellate_with_neuromaps", fake_parcellate)

    extracted = extract_scan_data(
        "lh.func.gii",
        atlas="dk",
        hemisphere="left",
        regions="all",
        source_space="fsaverage",
        input_rh="rh.func.gii",
    )

    np.testing.assert_allclose(extracted.values, expected)
    assert extracted.source_kind == "surface"
    assert seen["data"] == ("lh.func.gii", "rh.func.gii")
    assert seen["source_space"] == "fsaverage"
    assert seen["selection"].atlas.id == "dk"
    assert seen["selection"].n_regions == dk_left_selection.n_regions


def test_extract_scan_data_rejects_single_gifti_file():
    with pytest.raises(InputDataError, match="both left and right hemisphere files"):
        extract_scan_data("lh.func.gii", atlas="dk", hemisphere="left", regions="all")


def test_extract_scan_data_rejects_non_mni_volume_without_neuromaps(nifti_factory):
    path = nifti_factory(np.zeros((8, 8, 8), dtype=np.float32), name="native_space_scan.nii.gz")

    with pytest.raises(InputDataError, match="Non-MNI volumetric data require neuromaps-based resampling"):
        extract_scan_data(
            path,
            atlas="dk",
            hemisphere="left",
            regions="all",
            source_space="fsaverage",
            prefer_neuromaps=False,
        )


def test_extract_scan_data_uses_neuromaps_for_declared_non_mni_volume(monkeypatch, nifti_factory, dk_left_selection):
    path = nifti_factory(np.zeros((8, 8, 8), dtype=np.float32), name="cross_space_scan.nii.gz")
    expected = np.linspace(-1.0, 1.0, dk_left_selection.n_regions, dtype=float)
    seen = {}

    def fake_parcellate(data, selection, *, source_space, hemi=None):
        seen["data"] = data
        seen["selection"] = selection
        seen["source_space"] = source_space
        seen["hemi"] = hemi
        return expected

    monkeypatch.setattr(scan_module, "parcellate_with_neuromaps", fake_parcellate)

    extracted = extract_scan_data(
        path,
        atlas="dk",
        hemisphere="left",
        regions="all",
        source_space="fsaverage",
        prefer_neuromaps=True,
    )

    np.testing.assert_allclose(extracted.values, expected)
    assert extracted.source_kind == "volume"
    assert seen["data"] == path.as_posix()
    assert seen["selection"].atlas.id == "dk"
    assert seen["source_space"] == "fsaverage"


def test_region_means_ignores_nans_and_returns_nans_without_overlap():
    atlas_data = np.array(
        [
            [[1, 1], [2, 2]],
            [[0, 0], [0, 0]],
        ],
        dtype=np.int32,
    )
    image_data = np.array(
        [
            [[1.0, 2.0], [np.nan, 6.0]],
            [[5.0, 7.0], [9.0, 11.0]],
        ],
        dtype=float,
    )

    means = volume_inputs.region_means(image_data, atlas_data, np.array([1, 2, 3], dtype=int))
    no_overlap = volume_inputs.region_means(
        np.full((2, 2, 2), np.nan, dtype=float),
        atlas_data,
        np.array([1, 2], dtype=int),
    )

    np.testing.assert_allclose(means[:2], np.array([1.5, 6.0], dtype=float))
    assert np.isnan(means[2])
    assert np.isnan(no_overlap).all()


def test_preferred_resolution_falls_back_when_primary_grid_is_missing(dk_left_selection):
    atlas = replace(dk_left_selection.atlas, volume_1mm_path=None)
    fallback_selection = replace(dk_left_selection, atlas=atlas)
    image = nib.Nifti1Image(np.zeros((5, 5, 5), dtype=np.float32), np.diag([1.0, 1.0, 1.0, 1.0]))

    assert volume_inputs.preferred_resolution(image, fallback_selection) == "2mm"


def test_parcellate_with_neuromaps_uses_volume_branch_and_indexes_ids(
    monkeypatch,
    dk_left_selection,
    fake_parcellater_factory,
):
    atlas_ids = dk_left_selection.labels["id"].astype(int).to_numpy()
    fake_parcellater, calls = fake_parcellater_factory(np.arange(1, atlas_ids.max() + 4, dtype=float))
    monkeypatch.setattr(neuromaps_inputs, "import_neuromaps", lambda: fake_parcellater)

    values = neuromaps_inputs.parcellate_with_neuromaps(
        "scan.nii.gz",
        dk_left_selection,
        source_space="MNI152",
    )

    np.testing.assert_array_equal(values, atlas_ids.astype(float))
    assert calls["space"] == "MNI152"
    assert calls["resampling_target"] == "parcellation"
    assert calls["fit"] == {
        "data": "scan.nii.gz",
        "source_space": "MNI152",
        "ignore_background_data": True,
        "hemi": None,
    }


def test_parcellate_with_neuromaps_rejects_vectors_too_short_for_surface_ids(
    monkeypatch,
    dk_both_selection,
    fake_parcellater_factory,
):
    atlas_ids = dk_both_selection.labels["id"].astype(int).to_numpy()
    fake_parcellater, calls = fake_parcellater_factory(np.arange(1, atlas_ids.max(), dtype=float))
    monkeypatch.setattr(neuromaps_inputs, "import_neuromaps", lambda: fake_parcellater)
    monkeypatch.setattr(neuromaps_inputs, "load_surface_parcellation", lambda atlas, hemi: "parcel")

    with pytest.raises(InputDataError, match="too short for atlas ids up to"):
        neuromaps_inputs.parcellate_with_neuromaps(
            ("lh.func.gii", "rh.func.gii"),
            dk_both_selection,
            source_space="fsaverage",
        )
    assert calls["space"] == dk_both_selection.atlas.surface_space


def test_parcellate_nifti_direct_raises_alignment_error_without_resample(nifti_factory, dk_left_selection):
    path = nifti_factory(np.zeros((32, 40, 36), dtype=np.float32), name="native_like_T1w.nii.gz")

    with pytest.raises(InputDataError, match="Register the image to MNI152 first"):
        volume_inputs.parcellate_nifti_direct(path, dk_left_selection, allow_resample=False)
