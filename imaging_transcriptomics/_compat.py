from __future__ import annotations

import warnings
from contextlib import contextmanager


_PKG_RESOURCES_WARNING = (
    r"pkg_resources is deprecated as an API\. See "
    r"https://setuptools\.pypa\.io/en/latest/pkg_resources\.html\."
)


@contextmanager
def suppress_pkg_resources_deprecation():
    """Hide third-party pkg_resources deprecation warnings during optional imports."""

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=_PKG_RESOURCES_WARNING,
            category=UserWarning,
        )
        yield
