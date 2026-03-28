"""Public facade for result-bundle serialization helpers."""

from .outputs.metadata import metadata_dict, write_metadata_json
from .outputs.write import read_optional_table, run_to_tables, write_result_bundle

__all__ = [
    "metadata_dict",
    "read_optional_table",
    "run_to_tables",
    "write_metadata_json",
    "write_result_bundle",
]
