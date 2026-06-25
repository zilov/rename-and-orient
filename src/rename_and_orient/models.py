"""Data models for rename-and-orient."""
from dataclasses import dataclass


@dataclass
class PAFRecord:
    """Single PAF alignment record."""

    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    num_matches: int
    alignment_length: int
    mapping_quality: int

    @property
    def alignment_block_length(self) -> int:
        """Length of the alignment block on query."""
        return self.query_end - self.query_start


@dataclass
class ChromosomeMapping:
    """Mapping information for a single chromosome."""

    query_name: str
    query_length: int
    target_name: str
    total_alignment_length: int
    coverage: float
    plus_strand_length: int
    minus_strand_length: int
    needs_reverse_complement: bool
    target_prefix: str = "chr_"

    @property
    def target_suffix(self) -> str:
        """Extract suffix from target name (e.g., 'chr_5' -> '5', 'chrW' -> 'W')."""
        if self.target_name.startswith(self.target_prefix):
            return self.target_name[len(self.target_prefix):]
        return self.target_name


@dataclass
class UnlocMapping:
    """Mapping information for an unlocalized contig (SUPER_N_unloc_M)."""

    contig_name: str
    parent_chromosome: str
    unloc_number: int
    needs_reverse_complement: bool


@dataclass
class FinalChromosomeAssignment:
    """Final assignment for a chromosome after conflict resolution."""

    original_name: str
    new_name: str
    new_suffix: str
    needs_reverse_complement: bool
    is_sex_chromosome: bool
