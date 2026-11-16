"""rename-and-orient: rename and orient chromosomes based on PAF alignment."""
import logging

from .alignment import (
    build_chromosome_mappings,
    calculate_target_alignments,
    merge_intervals,
)
from .assignments import (
    build_unloc_mappings,
    resolve_chromosome_assignments,
    sort_assignments_for_output,
)
from .fasta import read_fasta, reverse_complement, write_fasta
from .mapping_table import load_mapping_table_assignments
from .models import (
    ChromosomeMapping,
    FinalChromosomeAssignment,
    PAFRecord,
    UnlocMapping,
)
from .names import (
    _HAP_SUFFIX_RE,
    SEX_CHROMOSOME_SUFFIXES,
    extract_chromosome_suffix,
    is_autosome_suffix,
    is_hap_contig,
    is_sex_chromosome_suffix,
    is_unloc_contig,
    parse_unloc_name,
    strip_hap_suffix,
)
from .paf import (
    detect_reference_prefix,
    filter_paf_records,
    parse_paf,
    validate_paf_fasta_consistency,
)
from .report import write_chromosome_list

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = [
    "PAFRecord", "ChromosomeMapping", "FinalChromosomeAssignment", "UnlocMapping",
    "SEX_CHROMOSOME_SUFFIXES", "_HAP_SUFFIX_RE",
    "strip_hap_suffix", "is_hap_contig", "is_unloc_contig",
    "is_sex_chromosome_suffix", "is_autosome_suffix",
    "extract_chromosome_suffix", "parse_unloc_name",
    "read_fasta", "parse_paf", "write_fasta", "write_chromosome_list", "reverse_complement",
    "build_chromosome_mappings", "filter_paf_records", "validate_paf_fasta_consistency",
    "detect_reference_prefix", "merge_intervals", "calculate_target_alignments",
    "resolve_chromosome_assignments", "sort_assignments_for_output", "build_unloc_mappings",
    "load_mapping_table_assignments",
]
