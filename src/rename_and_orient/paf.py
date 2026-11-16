"""PAF parsing, filtering and consistency checks."""
import logging
import re
from collections import defaultdict
from collections.abc import Collection
from pathlib import Path

from .models import PAFRecord

logger = logging.getLogger(__name__)


def parse_paf(paf_path: Path) -> list[PAFRecord]:
    """Parse PAF file.

    PAF format (tab-separated):
    1. Query sequence name
    2. Query sequence length
    3. Query start (0-based)
    4. Query end (0-based, open)
    5. Strand ('+' or '-')
    6. Target sequence name
    7. Target sequence length
    8. Target start (0-based)
    9. Target end (0-based, open)
    10. Number of matching bases
    11. Alignment block length
    12. Mapping quality (0-255)

    Lines with fewer than 12 tab-separated fields are skipped silently.

    Args:
        paf_path: Path to PAF file

    Returns:
        List of PAFRecord objects
    """
    records = []

    with open(paf_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')
            if len(fields) < 12:
                continue

            record = PAFRecord(
                query_name=fields[0],
                query_length=int(fields[1]),
                query_start=int(fields[2]),
                query_end=int(fields[3]),
                strand=fields[4],
                target_name=fields[5],
                target_length=int(fields[6]),
                target_start=int(fields[7]),
                target_end=int(fields[8]),
                num_matches=int(fields[9]),
                alignment_length=int(fields[10]),
                mapping_quality=int(fields[11])
            )
            records.append(record)

    return records


def detect_reference_prefix(records: list[PAFRecord]) -> str:
    """Auto-detect reference chromosome prefix from PAF target names.

    Sums query alignment length per prefix and returns the winner, which
    correctly favours main chromosomes (chr_, SUPER_) over many short-hit
    unplaced contigs (NW_, NT_). Logs the result to stderr.

    Args:
        records: PAF records whose target names are inspected.

    Returns:
        The detected reference chromosome prefix.

    Raises:
        ValueError: If no alphabetic prefix is found -- use
            --reference-chromosome-prefix to specify it explicitly.
    """
    prefix_aln_length: dict[str, int] = defaultdict(int)

    for record in records:
        m = re.match(r'^([A-Za-z_\.]+)', record.target_name)
        if m:
            prefix_aln_length[m.group(1)] += record.query_end - record.query_start

    if not prefix_aln_length:
        raise ValueError(
            "Could not detect reference chromosome prefix from PAF file: "
            "no target names with a recognisable alphabetic prefix were found. "
            "Please specify --reference-chromosome-prefix explicitly."
        )

    detected = max(prefix_aln_length, key=lambda p: prefix_aln_length[p])
    logger.info(
        "Auto-detected reference prefix: '%s' (total aligned bases: %s)",
        detected, f"{prefix_aln_length[detected]:,}",
    )
    return detected


def filter_paf_records(
    records: list[PAFRecord],
    query_chromosome_prefix: str = "SUPER_",
    target_prefix: str | None = None
) -> tuple[list[PAFRecord], str]:
    """Filter PAF records to keep only those matching query and target prefixes.

    Auto-detects target prefix if not specified.

    Args:
        records: List of PAF records
        query_chromosome_prefix: Required prefix for query names
        target_prefix: Required prefix for target names (auto-detected if None)

    Returns:
        Tuple of (filtered records, detected target prefix)
    """
    if target_prefix is None:
        target_prefix = detect_reference_prefix(records)

    filtered = []
    for record in records:
        if (record.query_name.startswith(query_chromosome_prefix) and
                record.target_name.startswith(target_prefix)):
            filtered.append(record)
    return filtered, target_prefix


def validate_paf_fasta_consistency(
    paf_records: list[PAFRecord],
    fasta_sequences: Collection[str],
    query_chromosome_prefix: str = "SUPER_"
) -> tuple[bool, list[str], list[str]]:
    """Validate consistency between PAF and FASTA files for chromosomes with given prefix.

    Checks that all chromosomes with query_chromosome_prefix in PAF exist in FASTA
    and vice versa.

    Args:
        paf_records: List of PAF records
        fasta_sequences: Sequence names in the FASTA (any collection, e.g. dict or FastaIndex)
        query_chromosome_prefix: Prefix for chromosome names

    Returns:
        Tuple of:
        - is_valid: True if all chromosomes match
        - in_paf_not_fasta: List of chromosomes in PAF but not in FASTA
        - in_fasta_not_paf: List of chromosomes in FASTA but not in PAF
    """
    paf_chromosomes = set()
    for record in paf_records:
        if (record.query_name.startswith(query_chromosome_prefix) and
                "_unloc_" not in record.query_name):
            paf_chromosomes.add(record.query_name)

    fasta_chromosomes = set()
    for name in fasta_sequences:
        if name.startswith(query_chromosome_prefix) and "_unloc_" not in name:
            fasta_chromosomes.add(name)

    in_paf_not_fasta = sorted(paf_chromosomes - fasta_chromosomes)
    in_fasta_not_paf = sorted(fasta_chromosomes - paf_chromosomes)

    is_valid = len(in_paf_not_fasta) == 0 and len(in_fasta_not_paf) == 0

    return is_valid, in_paf_not_fasta, in_fasta_not_paf
