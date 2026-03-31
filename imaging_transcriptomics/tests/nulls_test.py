from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from imaging_transcriptomics.exceptions import AtlasAssetError, NullModelError
from imaging_transcriptomics.models import ExtractedScan
import imaging_transcriptomics.nulls as spatial_nulls


def test_standardize_vector_handles_constant_and_nonconstant_inputs():
    constant = spatial_nulls._standardize_vector(np.array([3.0, 3.0, np.nan], dtype=float))
    varying = spatial_nulls._standardize_vector(np.array([1.0, 2.0, 3.0, np.nan], dtype=float))

    np.testing.assert_allclose(constant[:2], np.zeros(2, dtype=float))
    assert np.isnan(constant[2])
    assert abs(np.nanmean(varying)) < 1e-12
    assert np.isclose(np.nanstd(varying, ddof=1), 1.0)


def test_shuffle_within_groups_preserves_group_membership_and_singletons():
    values = np.array([1.0, 2.0, 10.0, 20.0, 30.0, 999.0], dtype=float)
    groups = np.array(["L", "L", "R", "R", "R", "B"], dtype=object)

    permuted = spatial_nulls.shuffle_within_groups(
        values,
        groups,
        n_permutations=8,
        rng=np.random.default_rng(42),
    )

    assert permuted.shape == (6, 8)
    for column in range(permuted.shape[1]):
        assert set(permuted[:2, column]) == {1.0, 2.0}
        assert set(permuted[2:5, column]) == {10.0, 20.0, 30.0}
        assert permuted[5, column] == 999.0


def test_normalize_null_maps_accepts_transposed_output_and_rejects_bad_shape():
    values_by_perm = np.arange(12, dtype=float).reshape(3, 4)

    np.testing.assert_array_equal(
        spatial_nulls._normalize_null_maps(values_by_perm, 3, 4, "vasa"),
        values_by_perm,
    )
    np.testing.assert_array_equal(
        spatial_nulls._normalize_null_maps(values_by_perm.T, 3, 4, "vasa"),
        values_by_perm,
    )
    with pytest.raises(NullModelError, match="unexpected shape"):
        spatial_nulls._normalize_null_maps(np.zeros((5, 5), dtype=float), 3, 4, "vasa")


def test_generate_surface_nulls_rejects_missing_surface_assets(dk_left_selection):
    atlas = replace(
        dk_left_selection.atlas,
        lh_surface_path=None,
        rh_surface_path=None,
        lh_annot_path=None,
        rh_annot_path=None,
    )
    broken_selection = replace(dk_left_selection, atlas=atlas)

    with pytest.raises(AtlasAssetError, match="does not have surface parcellation files"):
        spatial_nulls.generate_surface_nulls(
            np.linspace(-1.0, 1.0, 34, dtype=float),
            broken_selection,
            method="vasa",
            n_permutations=4,
            seed=1234,
        )


def test_generate_surface_nulls_normalizes_neuromaps_orientation(monkeypatch, dk_left_selection):
    calls = {}

    class FakeNulls:
        @staticmethod
        def vasa(**kwargs):
            calls.update(kwargs)
            n_perm = kwargs["n_perm"]
            n_values = kwargs["data"].shape[0]
            return np.arange(n_perm * n_values, dtype=float).reshape(n_perm, n_values)

    monkeypatch.setattr(spatial_nulls, "_import_neuromaps_nulls", lambda: FakeNulls)
    monkeypatch.setattr(spatial_nulls, "_surface_parcellation", lambda *args, **kwargs: "parcel")
    monkeypatch.setattr(spatial_nulls, "surface_geometry_paths", lambda *args, **kwargs: None)

    result = spatial_nulls.generate_surface_nulls(
        np.linspace(-1.0, 1.0, 34, dtype=float),
        dk_left_selection,
        method="vasa",
        n_permutations=4,
        seed=1234,
    )

    assert result.shape == (34, 4)
    np.testing.assert_array_equal(result[:, 0], np.arange(34, dtype=float))
    assert calls["atlas"] == dk_left_selection.atlas.surface_space
    assert calls["parcellation"] == "parcel"
    assert "surfaces" not in calls


def test_permute_scan_values_handles_non_cortical_only_inputs_without_surface_nulls(non_cortical_extracted):
    permuted, used_method = spatial_nulls.permute_scan_values(
        non_cortical_extracted,
        n_permutations=6,
        null_method="auto",
        seed=1234,
    )

    assert used_method == "auto"
    assert permuted.shape == (non_cortical_extracted.values.shape[0], 6)

    standardized = spatial_nulls._standardize_vector(non_cortical_extracted.values)
    groups = non_cortical_extracted.labels["hemisphere"].astype(str).to_numpy()
    for group in np.unique(groups):
        idx = np.flatnonzero(groups == group)
        for column in range(permuted.shape[1]):
            np.testing.assert_allclose(
                np.sort(permuted[idx, column]),
                np.sort(standardized[idx]),
            )


def test_permute_scan_values_random_is_seed_reproducible(dk_left_extracted):
    first, method_a = spatial_nulls.permute_scan_values(
        dk_left_extracted,
        n_permutations=5,
        null_method="random",
        seed=42,
    )
    second, method_b = spatial_nulls.permute_scan_values(
        dk_left_extracted,
        n_permutations=5,
        null_method="random",
        seed=42,
    )

    np.testing.assert_allclose(first, second)
    assert method_a == method_b == "random"


def test_permute_scan_values_uses_surface_nulls_for_cortex_and_group_shuffle_for_rest(monkeypatch, dk_both_selection):
    labels = dk_both_selection.labels.reset_index(drop=True)
    values = np.linspace(-2.0, 2.0, labels.shape[0], dtype=float)
    extracted = ExtractedScan(
        values=values,
        selection=dk_both_selection,
        source="array",
        source_kind="vector",
    )
    cortical_idx = np.flatnonzero(labels["structure"].to_numpy(dtype=object) == "cortex")
    non_cortical_idx = np.flatnonzero(labels["structure"].to_numpy(dtype=object) != "cortex")
    fake_cortical = np.full((cortical_idx.size, 3), 7.5, dtype=float)

    monkeypatch.setattr(spatial_nulls, "generate_surface_nulls", lambda *args, **kwargs: fake_cortical)

    permuted, used_method = spatial_nulls.permute_scan_values(
        extracted,
        n_permutations=3,
        null_method="vasa",
        seed=1234,
    )

    np.testing.assert_allclose(permuted[cortical_idx], fake_cortical)
    standardized = spatial_nulls._standardize_vector(values)
    groups = labels.iloc[non_cortical_idx]["hemisphere"].astype(str).to_numpy()
    for group in np.unique(groups):
        idx = non_cortical_idx[np.flatnonzero(groups == group)]
        for column in range(permuted.shape[1]):
            np.testing.assert_allclose(
                np.sort(permuted[idx, column]),
                np.sort(standardized[idx]),
            )
    assert used_method == "vasa"


def test_permute_scan_values_rejects_unknown_method(non_cortical_extracted):
    with pytest.raises(NullModelError, match="Unknown null method 'bogus'"):
        spatial_nulls.permute_scan_values(
            non_cortical_extracted,
            n_permutations=4,
            null_method="bogus",
        )
