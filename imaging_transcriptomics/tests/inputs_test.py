import numpy as np
import pytest
from pathlib import Path
import gseapy

from imaging_transcriptomics.inputs import (
    extract_average,
    load_gene_expression,
    load_gene_labels,
    read_scan,
    get_geneset
)
from imaging_transcriptomics.errors import InvalidSizeError, InvalidFormatError


# READ SCAN FUNCTION TESTS
def test_read_scan(tdata_dir):
    """Test the read scan function."""
    scan = read_scan(tdata_dir / "MNI152_T1_1mm.nii.gz")
    assert scan.shape == (182, 218, 182)
    assert scan.dtype == np.float64


def test_read_scan_errors(tdata_dir):
    """Test the errors given by the read scan function"""
    with pytest.raises(FileNotFoundError):
        read_scan(tdata_dir / "MNI152_T1_1mm.nii")
    with pytest.raises(InvalidFormatError):
        read_scan(tdata_dir / "wrong_format.txt")


# EXTRACT AVERAGE FUNCTION TESTS
def test_extract_average():
    """Test the extract average function."""
    scan = np.ones((182, 218, 182))
    average = extract_average(scan)
    assert average.dtype == np.float64
    assert average.shape == (41,)
    np.testing.assert_array_equal(average, np.ones(41))


def test_extract_average_errors():
    """Test the errors given by the extract average function."""
    with pytest.raises(InvalidSizeError):
        extract_average(np.ones((182, 218, 182, 1)))
    with pytest.raises(InvalidSizeError):
        extract_average(np.ones((91, 102, 91)))


# LOAD GENE EXPRESSION FUNCTION TESTS
def test_load_gene_expression():
    """Test the load gene expression function."""
    expression = load_gene_expression(regions="all")
    assert expression.shape == (41, 15633)
    assert expression.dtype == np.float64
    np.testing.assert_almost_equal(expression[0, 0], -.281, decimal=3)
    np.testing.assert_almost_equal(expression[34, 18], 1.199, decimal=3)
    np.testing.assert_almost_equal(expression[39, 13], -3.201, decimal=3)


def test_load_gene_expression_cort():
    expression = load_gene_expression(regions="cort")
    assert expression.shape == (34, 15633)
    assert expression.dtype == np.float64
    np.testing.assert_almost_equal(expression[0, 0], -.281, decimal=3)
    np.testing.assert_almost_equal(expression[32, 18], .473, decimal=3)


def test_load_gene_expression_errors():
    """Test the errors given by the load gene expression function."""
    with pytest.raises(ValueError):
        load_gene_expression(regions="wrong_region")


# LOAD GENE LABELS FUNCTION TESTS
def test_gene_labels_load():
    """Test that the genes labels are loaded correctly."""
    labels = load_gene_labels()
    assert labels.shape == (15633, 1)
    assert labels[78] == "ABHD6"
    assert labels[1635] == "C6orf106"
    assert "SLC7A10" in labels
    assert "audhd49b" not in labels
    assert "LOC102723968" in labels


# GENESET FUNCTION TESTS
def test_get_geneset():
    """Test the geneset function."""
    geneset = get_geneset("lake")
    assert isinstance(geneset, str)
    assert Path(geneset).exists()
    geneset = get_geneset("pooled")
    assert isinstance(geneset, str)
    assert Path(geneset).exists()
    geneset = get_geneset("POOLed")
    assert isinstance(geneset, str)
    assert Path(geneset).exists()
    geneset = get_geneset("GO_Biological_Process_2017")
    assert isinstance(geneset, str)
    assert geneset in gseapy.get_library_name()

