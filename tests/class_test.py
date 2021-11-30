import pytest
import numpy as np

from imaging_transcriptomics.base import Imt


def test_imt_init():
    imt = Imt(np.random.rand(41), method="corr")
    assert imt.method == "corr"
    assert imt.scan_data.shape == (41,)


def test_imt_init_error():
    with pytest.raises(ValueError):
        Imt(np.random.rand(41), method="foo")

    with pytest.raises(Exception):
        Imt(np.random.rand(41),
            method="corr",
            scan_data=np.random.rand(42),
            regions="cort+sub")
    with pytest.raises(Exception):
        Imt(np.random.rand(41),
            method="corr",
            scan_data=np.random.rand(35),
            regions="cort")


