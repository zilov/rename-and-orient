"""Deprecated shim: the contents of this module moved to dedicated modules.

* FASTA and sequences: :mod:`rename_and_orient.fasta`
* PAF parsing: :mod:`rename_and_orient.paf`
* Reports (summary, TSV, chromosome list): :mod:`rename_and_orient.report`
* Name sorting: :mod:`rename_and_orient.names`
"""
import warnings

from .fasta import (
    _RC_TABLE,
    calculate_genome_length,
    read_fasta,
    reverse_complement,
    write_fasta,
)
from .names import _unloc_output_name, natural_sort_key
from .paf import parse_paf
from .report import print_mapping_summary, save_mapping_tsv, write_chromosome_list

warnings.warn(
    "rename_and_orient.io is deprecated; import from rename_and_orient.fasta, "
    "rename_and_orient.paf, rename_and_orient.report or rename_and_orient.names instead",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "_RC_TABLE", "_unloc_output_name", "calculate_genome_length", "natural_sort_key",
    "parse_paf", "print_mapping_summary", "read_fasta", "reverse_complement",
    "save_mapping_tsv", "write_chromosome_list", "write_fasta",
]
