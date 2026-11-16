"""Chromosome assignment resolution and sorting.

Two similar-looking models are involved:

* ``ChromosomeMapping`` (from :mod:`rename_and_orient.alignment`) is the raw
  alignment result: which reference chromosome a query matches best, coverage,
  strand lengths and the orientation decision.
* ``FinalChromosomeAssignment`` (produced here) is the final naming decision:
  the new name/suffix a chromosome gets in the output, resolved from the
  mappings (conflicts, sex chromosomes, haplotype tags) and used for writing files.
"""
import logging
from collections import defaultdict
from collections.abc import Collection

from .models import ChromosomeMapping, FinalChromosomeAssignment, UnlocMapping
from .names import (
    _HAP_SUFFIX_RE,
    extract_chromosome_suffix,
    is_autosome_suffix,
    is_sex_chromosome_suffix,
    is_unloc_contig,
    parse_unloc_name,
    strip_hap_suffix,
)

logger = logging.getLogger(__name__)


def resolve_chromosome_assignments(
    mappings: list[ChromosomeMapping],
    sequences: Collection[str],
    query_chromosome_prefix: str = "SUPER_",
    output_prefix: str = "SUPER_"
) -> list[FinalChromosomeAssignment]:
    """Resolve chromosome assignments with conflict handling and sex chromosome logic.

    This function handles:
    1. Sex chromosome detection (W, Z, X, Y, Z1, Z2, etc.)
    2. Conflict resolution when multiple chromosomes map to same target
    3. Autosome -> sex chromosome mapping (reassign to max_autosome + N)
    4. Sex chromosome -> autosome mapping (skip the autosome number)
    5. Unmapped chromosomes keep original names (may cause mapped chr to shift)

    Args:
        mappings: List of ChromosomeMapping objects from PAF analysis
        sequences: Names of all sequences (any collection, e.g. dict or FastaIndex)
        query_chromosome_prefix: Input chromosome prefix (e.g., "SUPER_", "scaffold_")
        output_prefix: Output chromosome prefix (e.g., "SUPER_", "chr_", "chr", "")

    Returns:
        List of FinalChromosomeAssignment for all chromosomes
    """
    assignments = []

    all_chr_names = {name for name in sequences
                     if name.startswith(query_chromosome_prefix) and "_unloc_" not in name}
    mapped_names = {m.query_name for m in mappings if "_unloc_" not in m.query_name}
    unmapped_names = all_chr_names - mapped_names

    autosome_mappings = []
    sex_mappings = []

    for m in mappings:
        query_suffix = extract_chromosome_suffix(m.query_name, query_chromosome_prefix)
        if "_unloc_" in m.query_name:
            continue

        if is_sex_chromosome_suffix(query_suffix):
            sex_mappings.append(m)
        else:
            autosome_mappings.append(m)

    # Reserved numbers: numbers used by unmapped chromosomes.
    # These chromosomes keep their original names, so we can't assign these numbers.
    reserved_numbers = set()
    for name in unmapped_names:
        suffix = extract_chromosome_suffix(name, query_chromosome_prefix)
        if is_autosome_suffix(suffix):
            reserved_numbers.add(int(strip_hap_suffix(suffix)))
            logger.info(
                "Reserved: number %s (unmapped %s keeps original name)",
                strip_hap_suffix(suffix), name,
            )

    skipped_numbers = set()

    for m in sex_mappings:
        query_suffix = extract_chromosome_suffix(m.query_name, query_chromosome_prefix)
        target_suffix = m.target_suffix

        if is_autosome_suffix(target_suffix):
            skipped_numbers.add(int(target_suffix))
            logger.warning(
                "%s (sex chr) -> %s (autosome): number %s will be skipped",
                m.query_name, m.target_name, target_suffix,
            )

        assignment = FinalChromosomeAssignment(
            original_name=m.query_name,
            new_name=f"{output_prefix}{query_suffix}",
            new_suffix=strip_hap_suffix(query_suffix),
            needs_reverse_complement=m.needs_reverse_complement,
            is_sex_chromosome=True
        )
        assignments.append(assignment)

    target_to_autosomes = defaultdict(list)
    autosomes_to_sex_target = []

    for m in autosome_mappings:
        target_suffix = m.target_suffix

        if is_sex_chromosome_suffix(target_suffix):
            autosomes_to_sex_target.append(m)
            logger.warning(
                "%s (autosome) -> %s (sex chr): will be reassigned",
                m.query_name, m.target_name,
            )
        elif is_autosome_suffix(target_suffix):
            target_num = int(target_suffix)
            target_to_autosomes[target_num].append((m, m.total_alignment_length))

    assigned_numbers = set()
    autosome_assignments = []
    deferred_autosomes = []

    for target_num, candidates in target_to_autosomes.items():
        if target_num in reserved_numbers:
            deferred_autosomes.extend([m for m, _ in candidates])
            deferred = [m.query_name for m, _ in candidates]
            logger.debug("Number %s reserved - deferring: %s", target_num, deferred)
            continue

        if target_num in skipped_numbers:
            deferred_autosomes.extend([m for m, _ in candidates])
            continue

        candidates.sort(key=lambda x: x[1], reverse=True)

        winner, _ = candidates[0]
        autosome_assignments.append((winner, target_num))
        assigned_numbers.add(target_num)

        for m, _ in candidates[1:]:
            deferred_autosomes.append(m)
            logger.info(
                "Conflict: %s lost to %s for number %s", m.query_name, winner.query_name, target_num
            )

    all_used_numbers = assigned_numbers | reserved_numbers
    max_autosome = max(all_used_numbers) if all_used_numbers else 0

    next_available = max_autosome + 1
    unavailable = assigned_numbers | reserved_numbers | skipped_numbers

    for m in deferred_autosomes:
        while next_available in unavailable:
            next_available += 1
        autosome_assignments.append((m, next_available))
        assigned_numbers.add(next_available)
        unavailable.add(next_available)
        logger.warning("Reassigned: %s -> %s%s", m.query_name, output_prefix, next_available)
        next_available += 1

    for m in autosomes_to_sex_target:
        while next_available in unavailable:
            next_available += 1
        autosome_assignments.append((m, next_available))
        assigned_numbers.add(next_available)
        unavailable.add(next_available)
        logger.warning(
            "Reassigned (was sex target): %s -> %s%s", m.query_name, output_prefix, next_available
        )
        next_available += 1

    for m, num in autosome_assignments:
        hap_match = _HAP_SUFFIX_RE.search(m.query_name)
        hap_str = hap_match.group(0) if hap_match else ""
        assignment = FinalChromosomeAssignment(
            original_name=m.query_name,
            new_name=f"{output_prefix}{num}{hap_str}",
            new_suffix=str(num),
            needs_reverse_complement=m.needs_reverse_complement,
            is_sex_chromosome=False
        )
        assignments.append(assignment)

    for name in unmapped_names:
        suffix = extract_chromosome_suffix(name, query_chromosome_prefix)
        bare_suffix = strip_hap_suffix(suffix)
        assignment = FinalChromosomeAssignment(
            original_name=name,
            new_name=f"{output_prefix}{suffix}",
            new_suffix=bare_suffix,
            needs_reverse_complement=False,
            is_sex_chromosome=is_sex_chromosome_suffix(suffix)
        )
        assignments.append(assignment)
        logger.info("Unmapped: %s -> %s%s (keeping orientation)", name, output_prefix, suffix)

    return assignments


def sort_assignments_for_output(
    assignments: list[FinalChromosomeAssignment]
) -> list[FinalChromosomeAssignment]:
    """Sort assignments for output: autosomes by number, then sex chromosomes alphabetically.

    Args:
        assignments: List of FinalChromosomeAssignment

    Returns:
        Sorted list
    """
    autosomes = [a for a in assignments if not a.is_sex_chromosome]
    sex_chrs = [a for a in assignments if a.is_sex_chromosome]

    # non-numeric suffixes (1A, MT) sort after numeric ones instead of crashing
    autosomes.sort(
        key=lambda x: (0, int(x.new_suffix)) if x.new_suffix.isdigit() else (1, x.new_suffix)
    )
    sex_chrs.sort(key=lambda x: x.new_suffix)

    return autosomes + sex_chrs


def build_unloc_mappings(
    sequences: Collection[str],
    chromosome_mappings: list[ChromosomeMapping]
) -> list[UnlocMapping]:
    """Build mappings for unlocalized contigs based on parent chromosome mappings.

    Unloc contigs (SUPER_N_unloc_M) inherit orientation from their parent
    chromosome (SUPER_N). They don't need alignment-based mapping.

    Args:
        sequences: Names of all sequences (any collection, e.g. dict or FastaIndex)
        chromosome_mappings: List of chromosome mappings

    Returns:
        List of UnlocMapping objects
    """
    parent_rc_lookup = {m.query_name: m.needs_reverse_complement for m in chromosome_mappings}

    unloc_mappings = []

    for seq_name in sequences:
        if not is_unloc_contig(seq_name):
            continue

        parent_chr, unloc_num = parse_unloc_name(seq_name)

        needs_rc = parent_rc_lookup.get(parent_chr, False)

        unloc_mapping = UnlocMapping(
            contig_name=seq_name,
            parent_chromosome=parent_chr,
            unloc_number=unloc_num,
            needs_reverse_complement=needs_rc
        )
        unloc_mappings.append(unloc_mapping)

    return unloc_mappings
