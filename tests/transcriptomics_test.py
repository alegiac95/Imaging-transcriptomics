import numpy as np
import pytest
import imaging_transcriptomics as imt


def test_init_basic():
    test_param = {"n_components": 3, "variance": None}
    data = np.ones(41)
    t = imt.ImagingTranscriptomics(data, **test_param)
    t_z_scores = np.empty(41)
    t_z_scores[:] = np.NaN
    assert t.n_components == 3
    assert t.var is None
    np.testing.assert_equal(t.zscore_data, t_z_scores)

