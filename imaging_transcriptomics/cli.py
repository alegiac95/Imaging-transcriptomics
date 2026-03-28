#!/usr/bin/env python3

"""Public command-line entrypoint for the imaging transcriptomics toolbox."""

from __future__ import annotations

from ._logging import configure_cli_logging
from .cli_support.parser import build_parser, parse_cmdline
from .cli_support.runners import dispatch_command, resolve_run_gsea as _resolve_run_gsea

__all__ = ["build_parser", "parse_cmdline", "main", "_resolve_run_gsea"]


def main():
    """Run the command-line interface."""

    parsed = parse_cmdline()
    configure_cli_logging()
    dispatch_command(parsed)


if __name__ == "__main__":
    main()
