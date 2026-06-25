"""rename-and-orient: rename and orient chromosomes based on PAF alignment."""
from .alignment import (
    build_chromosome_mappings,
    calculate_target_alignments,
    detect_reference_prefix,
    filter_paf_records,
    merge_intervals,
    validate_paf_fasta_consistency,
)
from .assignments import (
    build_unloc_mappings,
    resolve_chromosome_assignments,
    sort_assignments_for_output,
)
from .io import (
    parse_paf,
    read_fasta,
    reverse_complement,
    write_chromosome_list,
    write_fasta,
)
from .mapping_table import load_mapping_table_assignments
from .models import (
    ChromosomeMapping,
    FinalChromosomeAssignment,
    PAFRecord,
    UnlocMapping,
)
from .names import (
    SEX_CHROMOSOME_SUFFIXES,
    _HAP_SUFFIX_RE,
    extract_chromosome_suffix,
    is_autosome_suffix,
    is_hap_contig,
    is_sex_chromosome_suffix,
    is_unloc_contig,
    parse_unloc_name,
    strip_hap_suffix,
)

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
