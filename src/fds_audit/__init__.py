"""fds-audit — static checker for silent geometry failures in FDS input files."""

from .domain_check import (
    Box,
    FdsModel,
    Grid,
    Issue,
    check_file,
    check_model,
    parse_fds,
    parse_records,
    parse_text,
    report,
)

__version__ = "0.1.0"

__all__ = [
    "Box",
    "FdsModel",
    "Grid",
    "Issue",
    "check_file",
    "check_model",
    "parse_fds",
    "parse_records",
    "parse_text",
    "report",
    "__version__",
]
