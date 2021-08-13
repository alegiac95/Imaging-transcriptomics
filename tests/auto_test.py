import pytest

import imaging_transcriptomics as imt


def test_version():
    """Test that the version is imported correctly."""
    imported = dir(imt)
    assert "__version__" in imported


def test_modules_import():
    """Test that the submodules are all imported correctly."""
    imported = dir(imt)
    assert "bootstrap" in imported
    assert "errors" in imported
    assert "genes" in imported
    assert "inputs" in imported
    assert "permutations" in imported
    assert "reporting" in imported
    assert "transcriptomics" in imported


def test_functions_import():
    """Test that the functions are imported correctly."""
    imported = dir(imt)
    assert "bootstrap_pls" in imported
    assert "bootstrap_genes" in imported
    assert "load_gene_expression" in imported
    assert "load_gene_labels" in imported
    assert "get_components" in imported
    assert "read_scan" in imported
    assert "extract_average" in imported


def test_classes_import():
    """Test that the classes are imported correctly."""
    imported = dir(imt)
    assert "ImagingTranscriptomics" in imported
    assert "GeneResults" in imported


def test_not_in_module():
    """Test that an error is raised when trying to import a non existing modeule"""
    with pytest.raises(ImportError):
        from imaging_transcriptomics import outputs
