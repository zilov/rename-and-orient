"""Alignment analysis: grouping, interval merging, coverage and best-target choice."""
import logging
from collections import defaultdict

from .models import ChromosomeMapping, PAFRecord
from .names import is_unloc_contig
from .orientation import needs_reverse_complement_by_correlation

logger = logging.getLogger(__name__)


def group_alignments_by_query(records: list[PAFRecord]) -> dict[str, list[PAFRecord]]:
    """Group PAF records by query chromosome name.

    Args:
        records: List of PAF records

    Returns:
        Dictionary mapping query names to lists of PAF records
    """
    groups = defaultdict(list)
    for record in records:
        groups[record.query_name].append(record)
    return dict(groups)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping intervals into non-overlapping ones.

    Intervals are half-open (0-based, end exclusive), as in PAF coordinates.

    Args:
        intervals: List of (start, end) tuples

    Returns:
        Sorted list of non-overlapping (start, end) tuples
    """
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def calculate_target_alignments(records: list[PAFRecord]) -> dict[str, dict[str, int]]:
    """Calculate total alignment lengths for each query-target pair.

    Overlapping query intervals on the same target are merged before summing
    to avoid double-counting repeated alignments to the same region
    (e.g. when a repetitive reference chromosome attracts many duplicate hits).
    Strand totals use the same de-duplicated intervals per strand.

    Args:
        records: List of PAF records for a single query

    Returns:
        Dictionary mapping target names to dict with 'total', 'plus', 'minus' lengths
    """
    target_intervals: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: {'+': [], '-': []}
    )

    for record in records:
        strand = record.strand
        interval = (record.query_start, record.query_end)
        target_intervals[record.target_name][strand].append(interval)

    target_stats: dict[str, dict[str, int]] = {}
    for target_name, strands in target_intervals.items():
        plus_len = sum(e - s for s, e in merge_intervals(strands['+'
                                                                   ]))
        minus_len = sum(e - s for s, e in merge_intervals(strands['-']))
        all_intervals = strands['+'] + strands['-']
        total_len = sum(e - s for s, e in merge_intervals(all_intervals))
        target_stats[target_name] = {
            'total': total_len,
            'plus': plus_len,
            'minus': minus_len,
        }

    return target_stats


def determine_best_target(
    target_stats: dict[str, dict[str, int]]
) -> tuple[str | None, dict[str, int] | None]:
    """Determine the best target chromosome based on total alignment length.

    Args:
        target_stats: Dictionary from calculate_target_alignments

    Returns:
        Tuple of (best_target_name, stats_dict), or (None, None) if no target
        has a positive total alignment length
    """
    best_target = None
    best_stats = None
    max_length = 0

    for target_name, stats in target_stats.items():
        if stats['total'] > max_length:
            max_length = stats['total']
            best_target = target_name
            best_stats = stats

    return best_target, best_stats


def build_chromosome_mappings(
    records: list[PAFRecord],
    min_coverage: float = 0.5,
    target_prefix: str = "chr_"
) -> list[ChromosomeMapping]:
    """Build chromosome mappings from PAF records.

    For each query chromosome:
    1. Calculate alignment lengths to each target
    2. Select best target (maximum alignment length)
    3. Check coverage threshold
    4. Determine orientation based on strand statistics

    Args:
        records: Filtered PAF records
        min_coverage: Minimum coverage threshold (0.0-1.0)
        target_prefix: Prefix used in reference chromosome names (chr_ or chr)

    Returns:
        List of ChromosomeMapping objects
    """
    query_lengths = {}
    for record in records:
        query_lengths[record.query_name] = record.query_length

    groups = group_alignments_by_query(records)

    mappings = []

    for query_name, query_records in groups.items():
        query_length = query_lengths[query_name]

        target_stats = calculate_target_alignments(query_records)

        best_target, best_stats = determine_best_target(target_stats)

        if best_target is None or best_stats is None:
            continue

        coverage = best_stats['total'] / query_length

        if coverage < min_coverage:
            if not is_unloc_contig(query_name):
                logger.warning(
                    "%s -> %s coverage %.2f%% below threshold",
                    query_name, best_target, coverage * 100,
                )
            continue

        # Determine orientation using Pearson correlation method.
        # Correlation between query_mid and target_mid identifies the dominant
        # syntenic strand while ignoring large structural inversions.
        best_target_records = [r for r in query_records if r.target_name == best_target]
        target_length = best_target_records[0].target_length if best_target_records else 0
        needs_rc = needs_reverse_complement_by_correlation(best_target_records, target_length)

        mapping = ChromosomeMapping(
            query_name=query_name,
            query_length=query_length,
            target_name=best_target,
            total_alignment_length=best_stats['total'],
            coverage=coverage,
            plus_strand_length=best_stats['plus'],
            minus_strand_length=best_stats['minus'],
            needs_reverse_complement=needs_rc,
            target_prefix=target_prefix
        )
        mappings.append(mapping)

    return mappings


# Backward-compatible re-exports of names that moved to other modules.
from .paf import (  # noqa: E402, F401
    detect_reference_prefix,
    filter_paf_records,
    validate_paf_fasta_consistency,
)
