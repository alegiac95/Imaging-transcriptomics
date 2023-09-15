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


def test_from_scan_init(tdata_dir):
    scan = tdata_dir / "MNI152_T1_1mm.nii.gz"
    imt_instance = imt.ImagingTranscriptomics.from_scan(scan, regions="all",
                                                        method="corr")
    assert isinstance(imt_instance, imt.ImagingTranscriptomics)
    with pytest.raises(ValueError):
        imt.ImagingTranscriptomics.from_scan(scan, regions="all",
                                             method="pca")
    with pytest.raises(ValueError):
        imt.ImagingTranscriptomics.from_scan(scan, regions="none",
                                             method="corr")


def test_from_scan_errors(tdata_dir):
    scan = tdata_dir / "MNI152_T1_1mm_brain.nii.gz"
    with pytest.raises(FileNotFoundError):
        imt.ImagingTranscriptomics.from_scan(scan, regions="all",
                                             method="corr")


def test_from_file_init(tdata_dir):
    file = tdata_dir / "test_input.txt"
    imt_instance = imt.ImagingTranscriptomics.from_file(file, regions="all",
                                                        method="corr")
    assert isinstance(imt_instance, imt.ImagingTranscriptomics)
    with pytest.raises(ValueError):
        imt.ImagingTranscriptomics.from_file(tdata_dir, regions="all",
                                             method="corr")
    with pytest.raises(FileNotFoundError):
        imt.ImagingTranscriptomics.from_file(tdata_dir / "new_file.txt",
                                             regions="none", method="corr")


def test_permute_data():
    """Test the permutations method."""
    data = np.random.rand(41)
    imt_instance = imt.ImagingTranscriptomics(data, regions="all",
                                              method="pls", n_components=1)
    assert imt_instance._permutations is None
    imt_instance.permute_data()
    assert imt_instance._permutations is not None
    assert imt_instance._permutations.shape == (41, 1000)
    assert imt_instance._permutations.dtype == np.float64


def test_permute_data_with_cort():
    """Test the permutations method."""
    data = np.random.rand(34)
    imt_instance = imt.ImagingTranscriptomics(data, regions="cort",
                                              method="pls", n_components=1)
    imt_instance.permute_data()
    assert imt_instance._permutations.shape == (34, 1000)
    assert imt_instance._permutations.dtype == np.float64


def test_make_out_dir(tmpdir):
    """Test the make_out_dir method."""
    out_dir = tmpdir / "Imt_test_pls"
    imt_instance = imt.ImagingTranscriptomics(np.random.rand(41),
                                              regions="all", method="pls",
                                              n_components=1)
    imt_instance._make_output_dir(tmpdir, "test")
    assert out_dir.exists()


def test_atlas_transcriptomics():
    imt_instance = imt.ImagingTranscriptomics(np.random.rand(41),
                                              atlas="DK")
    assert imt_instance.atlas == "DK"
    assert imt_instance.atlas != "Schaefer_100"
    imt_instance = imt.ImagingTranscriptomics(np.random.rand(50),
                                              atlas="Schaefer_100")
    assert imt_instance.atlas != "DK"
    assert imt_instance.atlas == "Schaefer_100"


def test_permute_data_schaefer_100():
    """Test the permutations method."""
    data = np.random.rand(50)
    imt_instance = imt.ImagingTranscriptomics(data,
                                              atlas="Schaefer_100",
                                              regions="all",
                                              method="corr")
    assert imt_instance._permutations is None
    imt_instance.permute_data()
    assert imt_instance._permutations is not None
    assert imt_instance._permutations.shape == (50, 1000)
    assert imt_instance._permutations.dtype == np.float64


def test_shaefer_error():
    """Test the correct number of inputs."""
    data = np.random.rand(41)
    with pytest.raises(ValueError):
        imt.ImagingTranscriptomics(data, atlas="Schaefer_100",
                                   regions="all", method="corr")

def test_run_schaefer():
    """Test the run method."""
    data = np.random.rand(50)
    imt_instance = imt.ImagingTranscriptomics(data,
                                              atlas="Schaefer_100",
                                              regions="all",
                                              method="corr")
    imt_instance.run(save_res = False)
    assert imt_instance._permutations is not None
    assert imt_instance._permutations.shape == (50, 1000)
    assert imt_instance._permutations.dtype == np.float64
    assert imt_instance._method == "corr"
    assert imt_instance.atlas == "Schaefer_100"
    assert imt_instance._regions == "all"
    assert imt_instance.scan_data is not None
    assert imt_instance._subcortical is None
    assert imt_instance.zscore_data is not None
