from pathlib import Path
import numpy as np
import pytest
import imaging_transcriptomics as imt
from scipy.stats import zscore


# INITIALIZATION TESTS
def test_init_transcriptomics_corr():
    """
    Test the initialization of a transcriptomics object.
    """
    # Test with all regions (cort + sub)
    data = np.random.rand(41)
    imt_instance = imt.ImagingTranscriptomics(data, regions="all",
                                              method="corr")
    assert isinstance(imt_instance, imt.ImagingTranscriptomics)
    assert imt_instance._regions == "all"
    assert imt_instance._method == "corr"
    assert imt_instance.scan_data.shape == (41,)
    assert imt_instance.scan_data.dtype == np.float64
    np.testing.assert_array_almost_equal(imt_instance._cortical, zscore(data,
                                                                 ddof=1)[:34])
    assert imt_instance.zscore_data.shape == (41,)
    assert imt_instance._permutations is None
    # Test with only cortical regions
    data = np.random.rand(34)
    imt_instance = imt.ImagingTranscriptomics(data, regions="cort",
                                              method="corr")
    assert isinstance(imt_instance, imt.ImagingTranscriptomics)
    assert imt_instance._regions == "cort"
    assert imt_instance._method == "corr"
    assert imt_instance.scan_data.shape == (34,)
    assert imt_instance.scan_data.dtype == np.float64
    assert imt_instance._subcortical is None
    np.testing.assert_array_almost_equal(imt_instance._cortical, zscore(
        data, ddof=1))
    assert imt_instance.zscore_data.shape == (34,)
    assert imt_instance._permutations is None


def test_init_transcriptomics_pls():
    """
    Test the initialization of a transcriptomics object.
    """
    # Test with a all regions (cort + sub)
    data = np.random.rand(41)
    imt_instance = imt.ImagingTranscriptomics(data, regions="all",
                                              method="pls", n_components=1)
    assert isinstance(imt_instance, imt.ImagingTranscriptomics)
    assert imt_instance._regions == "all"
    assert imt_instance._method == "pls"
    assert imt_instance.scan_data.shape == (41,)
    assert imt_instance.scan_data.dtype == np.float64
    np.testing.assert_array_almost_equal(imt_instance._cortical, zscore(
        data, ddof=1)[:34])
    assert imt_instance.zscore_data.shape == (41,)
    assert imt_instance._permutations is None
    # Test with only cortical regions
    data = np.random.rand(41)
    imt_instance = imt.ImagingTranscriptomics(data, regions="cort",
                                              method="pls", n_components=1)
    assert isinstance(imt_instance, imt.ImagingTranscriptomics)
    assert imt_instance._regions == "cort"
    assert imt_instance._method == "pls"
    assert imt_instance.scan_data.shape == (41,)
    assert imt_instance.scan_data.dtype == np.float64
    assert imt_instance._subcortical is None
    assert imt_instance.zscore_data.shape == (41,)
    assert imt_instance._permutations is None


def test_wrong_method():
    """
    Test the initialization of a transcriptomics object with a wrong method.
    """
    data = np.random.rand(41)
    with pytest.raises(ValueError):
        imt.ImagingTranscriptomics(data, regions="all", method="wrong")


def test_missing_pls_argument():
    """
    Test the initialization of a transcriptomics object with a wrong method.
    """
    data = np.random.rand(41)
    with pytest.raises(ValueError):
        imt.ImagingTranscriptomics(data, regions="all", method="pls")
