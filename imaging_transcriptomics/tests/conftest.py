from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from imaging_transcriptomics import select_atlas_data
from imaging_transcriptomics.models import ExtractedScan


@pytest.fixture(scope="session")
def tdata_dir() -> Path:
    """Return the packaged test-data directory."""

    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def dk_left_selection():
    """Return the default left-hemisphere DK atlas selection used in tests."""

    return select_atlas_data(atlas="dk", hemisphere="left", regions="all")


@pytest.fixture(scope="session")
def dk_both_selection():
    """Return the bilateral DK atlas selection used in tests."""

    return select_atlas_data(atlas="dk", hemisphere="both", regions="all")


@pytest.fixture
def dk_left_vector(dk_left_selection) -> np.ndarray:
    """Build a deterministic regional vector aligned to the DK left selection."""

    return np.linspace(-1.0, 1.0, dk_left_selection.n_regions, dtype=float)


@pytest.fixture
def dk_left_extracted(dk_left_selection, dk_left_vector) -> ExtractedScan:
    """Return an extracted regional vector for tests that bypass scan I/O."""

    return ExtractedScan(
        values=dk_left_vector,
        selection=dk_left_selection,
        source="array",
        source_kind="vector",
    )


@pytest.fixture
def non_cortical_extracted(dk_both_selection) -> ExtractedScan:
    """Return a subcortical-only extracted vector for null-model edge cases."""

    labels = dk_both_selection.labels.loc[dk_both_selection.labels["structure"] != "cortex"].reset_index(drop=True)
    subset = replace(
        dk_both_selection,
        labels=labels,
        expression=dk_both_selection.expression.iloc[: labels.shape[0]].reset_index(drop=True),
    )
    values = np.linspace(-1.0, 1.0, labels.shape[0], dtype=float)
    return ExtractedScan(values=values, selection=subset, source="array", source_kind="vector")


@pytest.fixture
def nifti_factory(tmp_path):
    """Create NIfTI files with configurable data, affine, and filename."""

    def factory(
        data: np.ndarray,
        *,
        name: str = "image.nii.gz",
        affine: np.ndarray | None = None,
    ) -> Path:
        image = nib.Nifti1Image(
            np.asarray(data, dtype=np.float32),
            np.eye(4, dtype=float) if affine is None else np.asarray(affine, dtype=float),
        )
        path = tmp_path / name
        nib.save(image, path)
        return path

    return factory


@pytest.fixture
def fake_parcellater_factory():
    """Build a fake neuromaps Parcellater class and capture constructor calls."""

    def factory(return_values) -> tuple[type, dict[str, object]]:
        calls: dict[str, object] = {}

        class FakeParcellater:
            def __init__(self, parcellation, space, resampling_target):
                calls["parcellation"] = parcellation
                calls["space"] = space
                calls["resampling_target"] = resampling_target

            def fit_transform(self, data, source_space, ignore_background_data=True, hemi=None):
                calls["fit"] = {
                    "data": data,
                    "source_space": source_space,
                    "ignore_background_data": ignore_background_data,
                    "hemi": hemi,
                }
                return np.asarray(return_values, dtype=float)

        return FakeParcellater, calls

    return factory
