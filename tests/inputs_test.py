import numpy as np
import pytest

from imaging_transcriptomics.inputs import extract_average, get_components, load_gene_expression, load_gene_labels


# def test_read_scan():
#   """Test that an imaging scan is red correctly."""
#   pass


def test_get_components():
    """Test that the returned number of components is correct."""
    test_variance = np.array([0.5, 0.1, 0.18, 0.22])
    assert get_components(0.50, test_variance) == 1
    assert get_components(0.60, test_variance) == 2
    assert get_components(0.61, test_variance) == 3
    assert get_components(1.00, test_variance) == 4


def test_get_components_error():
    """Test if the get_components() function raises the correct error."""
    test_variance = np.array([0, 0, 0])
    with pytest.raises(ValueError):
        get_components(1.1, test_variance)
    with pytest.raises(ValueError):
        get_components(-0.1, test_variance)


def test_extract_average():
    """Test that the extracted average from a scan is correct."""
    test_matrix = np.ones((182, 218, 182))

    np.testing.assert_equal(extract_average(test_matrix), np.ones(41))


def test_extract_average_len():
    """Test that the length of the returned object has the correct length (41)."""
    test_matrix = np.zeros((182, 218, 182))
    assert len(extract_average(test_matrix)) == 41


def test_gene_expression_load():
    """Test that the genes expression data are loaded correctly."""
    genes = load_gene_expression()
    assert genes.shape == (41, 15633)


def test_gene_labels_load():
    """Test that the genes labels are loaded correctly."""
    labels = load_gene_labels()
    assert labels.shape == (15633, 1)
    assert labels[78] == "ABHD6"
    assert labels[1635] == "C6orf106"
    assert "SLC7A10" in labels
    assert "audhd49b" not in labels  # just some randoms string to mimic a gene.

