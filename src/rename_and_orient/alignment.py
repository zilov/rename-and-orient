"""Alignment analysis: PAF filtering, coverage, orientation, plotting."""
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from .models import ChromosomeMapping, PAFRecord
from .names import is_unloc_contig


def _pearson_r(xs: List[float], ys: List[float]) -> float:
    """
    Pearson correlation coefficient between two lists of floats.
    Returns 0.0 if fewer than 2 points or zero variance in either variable.
    """
    n = len(xs)
    if n < 2:
        return 0.0
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx = sum((x - mx) ** 2 for x in xs) ** 0.5
    sy = sum((y - my) ** 2 for y in ys) ** 0.5
    if sx == 0 or sy == 0:
        return 0.0
    return num / (sx * sy)


def needs_reverse_complement_by_correlation(
    records: List[PAFRecord],
    target_length: int,
) -> bool:
    """
    Decide whether a query chromosome needs reverse-complement using Pearson
    correlation between query midpoints and target midpoints.

    For '+' blocks: r = corr(query_mid, target_mid).
    For '-' blocks: r = corr(query_mid, target_length - target_mid).

    The strand whose r is higher wins. This correctly ignores large structural
    inversions (their blocks are scattered -> low r) while syntenic blocks
    form a clean diagonal -> high r.

    Args:
        records:       PAFRecord objects for one query->best_target pair.
        target_length: Length of the target chromosome.

    Returns:
        True  if reverse complement is needed,
        False if the chromosome is already correctly oriented.
    """
    plus_recs = [r for r in records if r.strand == '+']
    minus_recs = [r for r in records if r.strand == '-']

    r_plus = 0.0
    if len(plus_recs) >= 2:
        xs = [(r.query_start + r.query_end) / 2 for r in plus_recs]
        ys = [(r.target_start + r.target_end) / 2 for r in plus_recs]
        r_plus = _pearson_r(xs, ys)

    r_minus = 0.0
    if len(minus_recs) >= 2:
        xs = [(r.query_start + r.query_end) / 2 for r in minus_recs]
        ys = [target_length - (r.target_start + r.target_end) / 2 for r in minus_recs]
        r_minus = _pearson_r(xs, ys)

    # If correlation is uninformative for both strands (< 2 blocks each),
    # fall back to naive length comparison.
    if r_plus == 0.0 and r_minus == 0.0:
        plus_len = sum(r.query_end - r.query_start for r in plus_recs)
        minus_len = sum(r.query_end - r.query_start for r in minus_recs)
        return minus_len > plus_len

    return r_minus > r_plus


def detect_reference_prefix(records: List[PAFRecord]) -> str:
    """
    Auto-detect reference chromosome prefix from PAF target names.

    Sums query alignment length per prefix and returns the winner, which
    correctly favours main chromosomes (chr_, SUPER_) over many short-hit
    unplaced contigs (NW_, NT_). Logs the result to stderr.

    Raises ValueError if no alphabetic prefix is found -- use
    --reference-chromosome-prefix to specify it explicitly.
    """
    prefix_aln_length: Dict[str, int] = defaultdict(int)

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
    print(f"[detect_reference_prefix] Auto-detected reference prefix: '{detected}' "
          f"(total aligned bases: {prefix_aln_length[detected]:,})", file=sys.stderr)
    return detected


def filter_paf_records(
    records: List[PAFRecord],
    query_chromosome_prefix: str = "SUPER_",
    target_prefix: str = None
) -> Tuple[List[PAFRecord], str]:
    """
    Filter PAF records to keep only those matching query and target prefixes.
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
    paf_records: List[PAFRecord],
    fasta_sequences: Dict[str, str],
    query_chromosome_prefix: str = "SUPER_"
) -> Tuple[bool, List[str], List[str]]:
    """
    Validate consistency between PAF and FASTA files for chromosomes with given prefix.

    Checks that all chromosomes with query_chromosome_prefix in PAF exist in FASTA
    and vice versa.

    Args:
        paf_records: List of PAF records
        fasta_sequences: Dictionary of FASTA sequences
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
    for name in fasta_sequences.keys():
        if name.startswith(query_chromosome_prefix) and "_unloc_" not in name:
            fasta_chromosomes.add(name)

    in_paf_not_fasta = sorted(paf_chromosomes - fasta_chromosomes)
    in_fasta_not_paf = sorted(fasta_chromosomes - paf_chromosomes)

    is_valid = len(in_paf_not_fasta) == 0 and len(in_fasta_not_paf) == 0

    return is_valid, in_paf_not_fasta, in_fasta_not_paf


def group_alignments_by_query(records: List[PAFRecord]) -> Dict[str, List[PAFRecord]]:
    """
    Group PAF records by query chromosome name.

    Args:
        records: List of PAF records

    Returns:
        Dictionary mapping query names to lists of PAF records
    """
    groups = defaultdict(list)
    for record in records:
        groups[record.query_name].append(record)
    return dict(groups)


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping intervals into non-overlapping ones.

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


def calculate_target_alignments(records: List[PAFRecord]) -> Dict[str, Dict[str, int]]:
    """
    Calculate total alignment lengths for each query-target pair.

    Overlapping query intervals on the same target are merged before summing
    to avoid double-counting repeated alignments to the same region
    (e.g. when a repetitive reference chromosome attracts many duplicate hits).
    Strand totals use the same de-duplicated intervals per strand.

    Args:
        records: List of PAF records for a single query

    Returns:
        Dictionary mapping target names to dict with 'total', 'plus', 'minus' lengths
    """
    target_intervals: Dict[str, Dict[str, List[Tuple[int, int]]]] = defaultdict(
        lambda: {'+': [], '-': []}
    )

    for record in records:
        strand = record.strand
        interval = (record.query_start, record.query_end)
        target_intervals[record.target_name][strand].append(interval)

    target_stats: Dict[str, Dict[str, int]] = {}
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
    target_stats: Dict[str, Dict[str, int]]
) -> Tuple[str, Dict[str, int]]:
    """
    Determine the best target chromosome based on total alignment length.

    Args:
        target_stats: Dictionary from calculate_target_alignments

    Returns:
        Tuple of (best_target_name, stats_dict)
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
    records: List[PAFRecord],
    min_coverage: float = 0.5,
    target_prefix: str = "chr_"
) -> List[ChromosomeMapping]:
    """
    Build chromosome mappings from PAF records.

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

        if best_target is None:
            continue

        coverage = best_stats['total'] / query_length

        if coverage < min_coverage:
            if not is_unloc_contig(query_name):
                print(f"  Warning: {query_name} -> {best_target} coverage {coverage:.2%} below threshold")
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


def plot_chromosome_alignments(
    records: List[PAFRecord],
    mappings: List[ChromosomeMapping],
    output_dir: Path,
) -> None:
    """
    Generate scatter plots of PAF alignment blocks for each mapped chromosome.

    Left panel:  raw alignment blocks (query_mid vs target_mid), coloured by strand.
    Right panel: same data after applying the orientation decision (target axis
                 flipped when RC is needed), so a correct decision shows a clean
                 diagonal from bottom-left to top-right.

    Args:
        records:    Filtered PAF records.
        mappings:   Chromosome mappings with orientation decisions.
        output_dir: Directory where PNG files are written.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "  Warning: matplotlib not installed -- skipping plots. "
            "Install with: pip install matplotlib",
            file=sys.stderr,
        )
        return

    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    blocks_index: Dict[str, List[PAFRecord]] = defaultdict(list)
    for r in records:
        blocks_index[f"{r.query_name}->{r.target_name}"].append(r)

    for mapping in mappings:
        if is_unloc_contig(mapping.query_name):
            continue
        key = f"{mapping.query_name}->{mapping.target_name}"
        recs = blocks_index.get(key, [])
        if not recs:
            continue

        plus_recs = [r for r in recs if r.strand == "+"]
        minus_recs = [r for r in recs if r.strand == "-"]
        target_len = max(r.target_end for r in recs)
        needs_rc = mapping.needs_reverse_complement

        fig, axes = plt.subplots(1, 2, figsize=(16, 7))

        ax = axes[0]
        if plus_recs:
            ax.scatter(
                [(r.target_start + r.target_end) / 2e6 for r in plus_recs],
                [(r.query_start + r.query_end) / 2e6 for r in plus_recs],
                s=1, alpha=0.4, color="steelblue", label="+ strand",
            )
        if minus_recs:
            ax.scatter(
                [(r.target_start + r.target_end) / 2e6 for r in minus_recs],
                [(r.query_start + r.query_end) / 2e6 for r in minus_recs],
                s=1, alpha=0.4, color="firebrick", label="- strand",
            )
        ax.set_xlabel(f"{mapping.target_name} position (Mb)")
        ax.set_ylabel(f"{mapping.query_name} position (Mb)")
        ax.set_title("All alignment blocks")
        ax.legend(markerscale=8, loc="upper left")

        ax2 = axes[1]
        if needs_rc:
            if plus_recs:
                ax2.scatter(
                    [(target_len - (r.target_start + r.target_end) / 2) / 1e6 for r in plus_recs],
                    [(r.query_start + r.query_end) / 2e6 for r in plus_recs],
                    s=1, alpha=0.4, color="steelblue", label="+ strand",
                )
            if minus_recs:
                ax2.scatter(
                    [(target_len - (r.target_start + r.target_end) / 2) / 1e6 for r in minus_recs],
                    [(r.query_start + r.query_end) / 2e6 for r in minus_recs],
                    s=1, alpha=0.4, color="firebrick", label="- strand",
                )
            ax2.set_title("After RC (target axis flipped)")
        else:
            if plus_recs:
                ax2.scatter(
                    [(r.target_start + r.target_end) / 2e6 for r in plus_recs],
                    [(r.query_start + r.query_end) / 2e6 for r in plus_recs],
                    s=1, alpha=0.4, color="steelblue", label="+ strand",
                )
            if minus_recs:
                ax2.scatter(
                    [(r.target_start + r.target_end) / 2e6 for r in minus_recs],
                    [(r.query_start + r.query_end) / 2e6 for r in minus_recs],
                    s=1, alpha=0.4, color="firebrick", label="- strand",
                )
            ax2.set_title("As-is (no RC applied)")

        ax2.set_xlabel(f"{mapping.target_name} position (Mb)")
        ax2.set_ylabel(f"{mapping.query_name} position (Mb)")
        ax2.legend(markerscale=8, loc="upper left")

        rc_str = "RC" if needs_rc else "no RC"
        cov_str = f"{mapping.coverage:.1%}"
        fig.suptitle(
            f"{mapping.query_name}  ->  {mapping.target_name}   |   "
            f"decision: {rc_str}   coverage: {cov_str}   "
            f"+: {mapping.plus_strand_length:,}  -: {mapping.minus_strand_length:,}",
            fontsize=11, fontweight="bold",
        )
        fig.tight_layout()
        out_path = plots_dir / f"{mapping.query_name}_vs_{mapping.target_name}.png"
        fig.savefig(out_path, dpi=120)
        plt.close(fig)

    print(f"  Plots saved to: {plots_dir}")
