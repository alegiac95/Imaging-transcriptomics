import pytest
import numpy as np
from imaging_transcriptomics.inputs import extract_average, read_scan, get_components


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
