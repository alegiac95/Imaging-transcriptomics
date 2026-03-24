import pytest

import imaging_transcriptomics as imt


def test_version():
    imported = dir(imt)
    assert "__version__" in imported


def test_v2_functions_import():
    imported = dir(imt)
    assert "run_corr" in imported
    assert "run_pls" in imported
    assert "run_analysis" in imported
    assert "build_run_config" in imported


def test_v2_types_import():
    imported = dir(imt)
    assert "RunConfig" in imported
    assert "CorrelationResult" in imported
    assert "PLSResult" in imported


def test_legacy_surface_removed():
    exported = set(imt.__all__)
    assert "ImagingTranscriptomics" not in exported
    assert "GeneResults" not in exported
    assert "read_scan" not in exported
    assert "extract_average" not in exported
    assert "inputs" not in exported


def test_not_in_module():
    with pytest.raises(ImportError):
        from imaging_transcriptomics import outputs
