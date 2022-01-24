from pathlib import Path


# Custom Errors to throw at the user
class InvalidFormatError(Exception):
    """Exception raised when the format of one of the input files in not correct.

    Attributes:
        errorFile -- file that is not in the correct format.
        message -- optional user overridden error message to display.
    """
    def __init__(
        self,
        error_file,
        message="The provided file has an invalid format. Please use files "
                "in the .nii, .nii.gz format.",
    ):
        self.error_file = error_file
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}" \
               f" The error was caused by the file {self.error_file}."


class InvalidSizeError(Exception):
    """Exception raised when the size of the images is not correct.

    Attributes:
        * errorFile -- file with the wrong size
        * size -- size of the input image
        * message -- optional user defined error message
    """
    def __init__(
        self, shape=(182, 218, 182),
            message="The provided file has a wrong shape."
    ):
        self.message = message
        self.shape = shape
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} The file has shape: {self.shape}"


# Checks Decorators
class CheckPath:
    """Decorator to check if a path exists.

    In order to run the function decorated the path provided to the function
    has to exists, otherwise an error is raised.

    :raises FileNotFoundError:
    """
    def __init__(self, function):
        self.function = function

    def __call__(self, path, *args, **kwargs):
        if not Path(path).absolute().exists():
            raise FileNotFoundError
        return self.function(path, *args, **kwargs)


class CheckExtension:
    """Decorator to check the file extension of the input scan.

    Extension of the imaging scan has to be in NIfTI format (compressed or not)
     in order to run the function.

    :raises InvalidFormatError:
    """
    def __init__(self, function):
        self.function = function

    def __call__(self, path, *args, **kwargs):
        imaging_path = Path(path)
        if str().join(imaging_path.suffixes) not in [".nii", ".nii.gz"]:
            raise InvalidFormatError(path)
        return self.function(path, *args, **kwargs)


class CheckShape:
    """Decorator to check the correct matrix shape of the imaging scan.

    Shape of the matrix has to be 182x218x182 in order to run the function,
    otherwise raises and error.

    :raises InvalidSizeError:
    """
    def __init__(self, function):
        self.function = function

    def __call__(self, image, *args, **kwargs):
        if not image.shape == (182, 218, 182):
            raise InvalidSizeError(image.shape)
        return self.function(image, *args, **kwargs)


class CheckVariance:
    """Decorator to check that the variance is in the correct range of values.

    Target variance has to be in the range 0.0-1.0 (equivalent to 0-100%) in
    order to run the function.

    :raises ValueError:
    """
    def __init__(self, function):
        self.function = function

    def __call__(self, target_var, *args, **kwargs):
        if target_var < 0.0 or target_var > 1.0:
            raise ValueError
        return self.function(target_var, *args, **kwargs)
