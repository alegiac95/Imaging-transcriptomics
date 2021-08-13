import numpy as np
import pytest
import imaging_transcriptomics as imt


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
