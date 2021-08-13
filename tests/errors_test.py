from pathlib import Path
import pytest
import numpy as np

from imaging_transcriptomics.errors import (
    CheckExtension,
    InvalidFormatError,
    CheckShape,
    InvalidSizeError,
    CheckVariance,
    CheckPath
)


def test_check_extension():
    """Test that the CheckExtension decorator works"""
    @CheckExtension
    def function(path):
        """Dummy function to test the decorator."""
        return path
    test_path = Path(__file__).parent / "data"
    p = function(test_path / "anatomical.nii")
    assert str(p) == str(Path(__file__).parent / "data" / "anatomical.nii")


def test_check_extension_error():
    """Test that the CheckExtension decorator throws the correct error."""
    @CheckExtension
    def function(path):
        """Dummy function to test the decorator."""
        return path

    test_path = Path(__file__).parent / "data"
    with pytest.raises(InvalidFormatError) as ex:
        function(test_path / "wrong_format.txt")

    assert str(ex.value) == f"The provided file has an invalid format. Please use files in the .nii, .nii.gz format. " \
                            f"The error was caused by the file {test_path.absolute()}/wrong_format.txt."


def test_check_shape():
    """Test the CheckShape decorator."""
    matrix = np.zeros((182, 218, 182))

    @CheckShape
    def function(in_matrix):
        """Dummy function to test the decorator."""
        return in_matrix.shape
    assert function(matrix) == (182, 218, 182)


def test_check_shape_error():
    """Test the CheckShape decorator throws the correct error."""
    matrix = np.zeros((171, 230, 167))

    @CheckShape
    def function(in_matrix):
        """Dummy function to test the decorator."""
        return in_matrix.shape

    with pytest.raises(InvalidSizeError) as ex:
        function(matrix)

    assert str(ex.value) == "The provided file has a wrong shape. The file has shape: (171, 230, 167)"


def test_check_variance():
    """Test the CheckVariance decorator."""
    @CheckVariance
    def function(var):
        """Dummy function to test the decorator."""
        return var + 1

    assert function(0.1) == 1.1


def test_check_variance_error():
    """Test that the CheckVariance decorator throws the correct error."""
    @CheckVariance
    def function(var):
        """Dummy function to test the decorator."""
        return var + 1

    with pytest.raises(ValueError):
        function(1.2)
    with pytest.raises(ValueError):
        function(-0.1)


def test_check_path(tmp_path):
    @CheckPath
    def function(path):
        return str(path)

    assert function(tmp_path) == str(tmp_path)
    with pytest.raises(FileNotFoundError):
        function(tmp_path / "foo_bar")
