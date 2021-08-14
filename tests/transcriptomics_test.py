from pathlib import Path
import numpy as np
import pytest
import imaging_transcriptomics as imt


# INITIALIZATION TESTS
def test_init_basic():
    """Test basic initialisation method."""
    test_param = {"n_components": 3, "variance": None}
    data = np.ones(41)
    t = imt.ImagingTranscriptomics(data, **test_param)
    t_z_scores = np.empty(41)
    t_z_scores[:] = np.NaN
    assert t.n_components == 3
    assert t.var is None
    np.testing.assert_equal(t.zscore_data, t_z_scores)


def test_init_error_components():
    """Test that the initialisation with components in the wrong range raises a ValueError."""
    t_param_err = [16, 0]
    for test in t_param_err:
        with pytest.raises(ValueError):
            imt.ImagingTranscriptomics(np.ones(41), n_components=test)


def test_init_error_variance():
    """Test that variance in the wrong range or with wrong type raises a ValueError."""
    t_param_err = [-1.0, 101]
    for test in t_param_err:
        with pytest.raises(ValueError):
            imt.ImagingTranscriptomics(np.ones(41), variance=test)
    test = 'a'
    with pytest.raises(TypeError):
        imt.ImagingTranscriptomics(np.ones(41), variance=test)


def test_init_missing_attribute_error():
    """Test that an AttributeError is raised when both number of components and variance are missing."""
    with pytest.raises(AttributeError):
        imt.ImagingTranscriptomics(np.ones(41))


def test_init_length_error():
    """Test an AttributeError is raised when the length of input data is different from 41"""
    test_param = {"n_components": 3, "variance": None}
    with pytest.raises(AttributeError):
        imt.ImagingTranscriptomics(np.ones(43), **test_param)
        imt.ImagingTranscriptomics(np.ones(40), **test_param)


# METHODS TESTS
def test_permutations():
    """Test the permutations are computed correctly."""
    t_data = np.array([2.49176123, 2.12076098, 1.88912675, 2.29363057, 2.17108429,
                       2.44944779, 1.9532944, 2.13951822, 1.81947959, 1.58996705,
                       1.91860982, 2.30857561, 2.39706742, 2.03412347, 2.0920649,
                       1.89473161, 2.05717326, 1.20646305, 1.72044527, 2.0083166,
                       1.66318842, 2.06091217, 1.72413881, 2.33628019, 2.61411213,
                       1.807411, 1.96163793, 1.85169722, 2.11455623, 1.92936416,
                       1.28974378, 1.81579151, 2.66449885, 2.67599858, 1.13808303,
                       1.40784474, 2.70367057, 2.00515875, 2.49107748, 1.75756543,
                       2.29094877])
    test = imt.ImagingTranscriptomics(t_data, n_components=1)
    assert test.permuted is None
    test.permute_data(10)
    assert test.permuted.shape == (41, 10)
    assert np.unique(test.permuted, axis=0).shape[1] == 10


def test_saving_error():
    """Test that the save_permutations() method raises an error if no permutations are available."""
    test = imt.ImagingTranscriptomics(np.zeros(41), n_components=1)
    with pytest.raises(AttributeError):
        test.save_permutations(Path().cwd())
