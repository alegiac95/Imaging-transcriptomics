from __future__ import annotations

import logging
import logging.config
from pathlib import Path

import yaml


_LOGGING_CONFIGURED = False


def get_logger(name: str) -> logging.Logger:
    global _LOGGING_CONFIGURED
    if not _LOGGING_CONFIGURED:
        cfg_file_path = Path(__file__).parent / "log_config.yaml"
        with open(cfg_file_path, "r") as config_file:
            log_cfg = yaml.safe_load(config_file.read())
        logging.config.dictConfig(log_cfg)
        _LOGGING_CONFIGURED = True
    return logging.getLogger(name)
