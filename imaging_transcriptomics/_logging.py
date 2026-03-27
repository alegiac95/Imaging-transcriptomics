from __future__ import annotations

import logging

PACKAGE_LOGGER_NAME = "imaging_transcriptomics"
_PACKAGE_LOGGER = logging.getLogger(PACKAGE_LOGGER_NAME)
if not any(isinstance(handler, logging.NullHandler) for handler in _PACKAGE_LOGGER.handlers):
    _PACKAGE_LOGGER.addHandler(logging.NullHandler())

_CLI_LOGGING_CONFIGURED = False


def get_logger(name: str | None = None) -> logging.Logger:
    """Return a package logger without configuring global logging side effects."""

    if name is None:
        return _PACKAGE_LOGGER
    if name.startswith(PACKAGE_LOGGER_NAME):
        return logging.getLogger(name)
    return logging.getLogger(f"{PACKAGE_LOGGER_NAME}.{name}")


def configure_cli_logging(level: int = logging.INFO) -> logging.Logger:
    """Configure console logging for the command-line interface only."""

    global _CLI_LOGGING_CONFIGURED
    package_logger = _PACKAGE_LOGGER
    if not _CLI_LOGGING_CONFIGURED:
        formatter = logging.Formatter("%(asctime)s: %(message)s", "%Y-%m-%d %H:%M:%S")
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        package_logger.handlers = [handler]
        package_logger.propagate = False
        _CLI_LOGGING_CONFIGURED = True
    package_logger.setLevel(level)
    for handler in package_logger.handlers:
        handler.setLevel(level)
    return package_logger
