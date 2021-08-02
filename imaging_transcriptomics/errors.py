from pathlib import Path


# Custom Errors to throw at the user


class InvalidFormatError(Exception):
    """Exception raised when the format of one of the input files in not correct.

    Attributes:
        errorFile -- file that is not in the correct format.
        message -- optional user overridden error message to display.
    """

    def __init__(self, error_file, message="The provided file has an invalid format. Please use files in the .nii, "
                                           ".nii.gz format."):
        self.error_file = error_file
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} The error was caused by the file {self.error_file}."


class InvalidSizeError(Exception):
    """Exception raised when the size of the images is not correct.

    Attributes:
        * errorFile -- file with the wrong size
        * size -- size of the input image
        * message -- optional user defined error message
    """

    def __init__(self, error_file, size=(182, 218, 182), message="The provided file has a wrong size."):
        self.error_file = error_file
        self.message = message
        self.size = size
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} The file {self.error_file} has size: {self.size}"


# Checks Decorators
# TODO: Change decorators from functions to classes.


def check_path_exists(func):
    """Decorator function to check path exists.
    """
    def wrapper(path, *args, **kwargs):
        if Path(path).exists():
            func(path, *args, **kwargs)
        else:
            raise FileNotFoundError
    return wrapper


def check_correct_shape(func):
    """Decorator to check the correct matrix size of an imaging file.
    """
    def wrapper(image, *args, **kwargs):
        if image.shape == (182, 218, 182):
            func(image, *args, **kwargs)
        else:
            raise InvalidSizeError
    return wrapper


def check_extensions(func):
    """Decorator to check for the correct file extension of the imaging file.
    """
    def wrapper(path, *args, **kwargs):
        imaging_path = Path(path)
        correct_suffixes = [".nii", ".nii.gz"]
        if str().join(imaging_path.suffixes) not in correct_suffixes:
            raise InvalidFormatError
        else:
            func(path, *args, **kwargs)
    return wrapper


def check_var_in_range(func):
    """Decorator to check that the variance is in the right range"""
    def wrapper(target_var, *args, **kwargs):
        if 0.0 < target_var <= 1.0:
            func(target_var, *args, **kwargs)
        else:
            raise AssertionError
    return wrapper
