"""Chromosome assignment resolution and sorting."""
from collections import defaultdict
from typing import Dict, List, Tuple

from .models import ChromosomeMapping, FinalChromosomeAssignment, UnlocMapping
from .names import (
    _HAP_SUFFIX_RE,
    extract_autosome_number,
    extract_chromosome_suffix,
    is_autosome_suffix,
    is_sex_chromosome_suffix,
    is_unloc_contig,
    parse_unloc_name,
    strip_hap_suffix,
)


def resolve_chromosome_assignments(
    mappings: List[ChromosomeMapping],
    sequences: Dict[str, str],
    query_chromosome_prefix: str = "SUPER_",
    output_prefix: str = "SUPER_"
) -> Tuple[List[FinalChromosomeAssignment], Dict[str, bool]]:
    """
    Resolve chromosome assignments with conflict handling and sex chromosome logic.

    This function handles:
    1. Sex chromosome detection (W, Z, X, Y, Z1, Z2, etc.)
    2. Conflict resolution when multiple chromosomes map to same target
    3. Autosome -> sex chromosome mapping (reassign to max_autosome + N)
    4. Sex chromosome -> autosome mapping (skip the autosome number)
    5. Unmapped chromosomes keep original names (may cause mapped chr to shift)

    Args:
        mappings: List of ChromosomeMapping objects from PAF analysis
        sequences: Dictionary of all sequences (for finding unmapped ones)
        query_chromosome_prefix: Input chromosome prefix (e.g., "SUPER_", "scaffold_")
        output_prefix: Output chromosome prefix (e.g., "SUPER_", "chr_", "chr", "")

    Returns:
        Tuple of:
        - List of FinalChromosomeAssignment for all chromosomes
        - Dictionary mapping original names to needs_reverse_complement flag
    """
    assignments = []
    rc_lookup = {}

    all_chr_names = {name for name in sequences.keys()
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
            reserved_numbers.add(extract_autosome_number(suffix))
            print(f"  Reserved: number {extract_autosome_number(suffix)} (unmapped {name} keeps original name)")

    skipped_numbers = set()

    for m in sex_mappings:
        query_suffix = extract_chromosome_suffix(m.query_name, query_chromosome_prefix)
        target_suffix = m.target_suffix

        if is_autosome_suffix(target_suffix):
            skipped_numbers.add(extract_autosome_number(target_suffix))
            print(f"  Note: {m.query_name} (sex chr) -> {m.target_name} (autosome): "
                  f"number {extract_autosome_number(target_suffix)} will be skipped")

        assignment = FinalChromosomeAssignment(
            original_name=m.query_name,
            new_name=f"{output_prefix}{query_suffix}",
            new_suffix=strip_hap_suffix(query_suffix),
            needs_reverse_complement=m.needs_reverse_complement,
            is_sex_chromosome=True
        )
        assignments.append(assignment)
        rc_lookup[m.query_name] = m.needs_reverse_complement

    target_to_autosomes = defaultdict(list)
    autosomes_to_sex_target = []

    for m in autosome_mappings:
        target_suffix = m.target_suffix

        if is_sex_chromosome_suffix(target_suffix):
            autosomes_to_sex_target.append(m)
            print(f"  Note: {m.query_name} (autosome) -> {m.target_name} (sex chr): "
                  f"will be reassigned")
        elif is_autosome_suffix(target_suffix):
            # Key by full suffix string (e.g. '1A', '1B', '7') to handle subgenome letters
            target_to_autosomes[target_suffix].append((m, m.total_alignment_length))

    assigned_numbers = set()
    autosome_assignments = []
    deferred_autosomes = []

    for target_suffix_key, candidates in target_to_autosomes.items():
        target_num = extract_autosome_number(target_suffix_key)
        if target_num in reserved_numbers:
            deferred_autosomes.extend([m for m, _ in candidates])
            print(f"  Number {target_num} reserved - deferring: {[m.query_name for m, _ in candidates]}")
            continue

        if target_num in skipped_numbers:
            deferred_autosomes.extend([m for m, _ in candidates])
            continue

        candidates.sort(key=lambda x: x[1], reverse=True)

        winner, _ = candidates[0]
        autosome_assignments.append((winner, target_suffix_key))
        assigned_numbers.add(target_num)

        for m, _ in candidates[1:]:
            deferred_autosomes.append(m)
            print(f"  Conflict: {m.query_name} lost to {winner.query_name} for {target_suffix_key}")

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
        print(f"  Reassigned: {m.query_name} -> {output_prefix}{next_available}")
        next_available += 1

    for m in autosomes_to_sex_target:
        while next_available in unavailable:
            next_available += 1
        autosome_assignments.append((m, next_available))
        assigned_numbers.add(next_available)
        unavailable.add(next_available)
        print(f"  Reassigned (was sex target): {m.query_name} -> {output_prefix}{next_available}")
        next_available += 1

    for m, suffix_or_num in autosome_assignments:
        hap_match = _HAP_SUFFIX_RE.search(m.query_name)
        hap_str = hap_match.group(0) if hap_match else ""
        new_suffix = str(suffix_or_num)
        assignment = FinalChromosomeAssignment(
            original_name=m.query_name,
            new_name=f"{output_prefix}{new_suffix}{hap_str}",
            new_suffix=new_suffix,
            needs_reverse_complement=m.needs_reverse_complement,
            is_sex_chromosome=False
        )
        assignments.append(assignment)
        rc_lookup[m.query_name] = m.needs_reverse_complement

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
        rc_lookup[name] = False
        print(f"  Unmapped: {name} -> {output_prefix}{suffix} (keeping orientation)")

    return assignments, rc_lookup


def _autosome_sort_key(suffix: str) -> tuple:
    """Sort key for autosome suffixes: (numeric_part, subgenome_letter).
    Handles both plain numbers ('7') and subgenome-annotated ('7B', '1A').
    """
    import re as _re
    m = _re.match(r'^(\d+)([A-Za-z]*)$', suffix)
    if m:
        return (int(m.group(1)), m.group(2))
    try:
        return (int(suffix), '')
    except ValueError:
        return (0, suffix)


def sort_assignments_for_output(
    assignments: List[FinalChromosomeAssignment]
) -> List[FinalChromosomeAssignment]:
    """
    Sort assignments for output: autosomes by number (then subgenome letter), then sex chromosomes alphabetically.

    Args:
        assignments: List of FinalChromosomeAssignment

    Returns:
        Sorted list
    """
    autosomes = [a for a in assignments if not a.is_sex_chromosome]
    sex_chrs = [a for a in assignments if a.is_sex_chromosome]

    autosomes.sort(key=lambda x: _autosome_sort_key(x.new_suffix))
    sex_chrs.sort(key=lambda x: x.new_suffix)

    return autosomes + sex_chrs


def build_unloc_mappings(
    sequences: Dict[str, str],
    chromosome_mappings: List[ChromosomeMapping]
) -> List[UnlocMapping]:
    """
    Build mappings for unlocalized contigs based on parent chromosome mappings.

    Unloc contigs (SUPER_N_unloc_M) inherit orientation from their parent
    chromosome (SUPER_N). They don't need alignment-based mapping.

    Args:
        sequences: Dictionary of all sequences (name -> sequence)
        chromosome_mappings: List of chromosome mappings

    Returns:
        List of UnlocMapping objects
    """
    parent_rc_lookup = {m.query_name: m.needs_reverse_complement for m in chromosome_mappings}

    unloc_mappings = []

    for seq_name in sequences.keys():
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
