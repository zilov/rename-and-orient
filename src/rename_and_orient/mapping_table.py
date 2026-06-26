"""Load pre-built mapping table for second-haplotype renaming."""
import re
from pathlib import Path
from typing import Dict, List, Tuple

from .models import FinalChromosomeAssignment, UnlocMapping
from .names import (
    _HAP_SUFFIX_RE,
    extract_chromosome_suffix,
    is_sex_chromosome_suffix,
    is_unloc_contig,
    parse_unloc_name,
    strip_hap_suffix,
)


def _suffix_from_renamed(renamed_to: str) -> str:
    """Extract the chromosomal suffix from a renamed_to value.

    Strips any _HAPJ suffix first, then strips the leading alphabetic/underscore
    prefix, returning just the chromosomal identifier (e.g. '3B', '1A', 'W', '7').
    Works regardless of which output prefix (SUPER_, chr_, chr, ...) was used.
    """
    base = strip_hap_suffix(renamed_to)
    m = re.match(r'^[A-Za-z_.]+(.+)$', base)
    return m.group(1) if m else base


def load_mapping_table_assignments(
    mapping_table_path: Path,
    sequences: Dict[str, str],
    query_chromosome_prefix: str = "SUPER_",
    output_prefix: str = "SUPER_",
) -> Tuple[List[FinalChromosomeAssignment], List[UnlocMapping]]:
    """
    Build assignments and unloc mappings from a pre-built mapping TSV
    (produced by a previous run on haplotype 1). Columns used: query,
    renamed_to, needs_reverse_complement. Unlocs present only in this
    haplotype are mapped via their parent chromosome; they are never RC'd.

    HAP suffix mismatch is handled transparently: if the table was built
    from _HAP1 sequences but the current FASTA contains _HAP2 sequences,
    the lookup is done on the base name (HAP suffix stripped) and the
    correct HAP suffix is re-applied to the output name.
    """
    table_map: Dict[str, Dict] = {}
    with open(mapping_table_path) as fh:
        header = None
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                missing = {"query", "renamed_to", "needs_reverse_complement"} - set(header)
                if missing:
                    raise ValueError(f"Mapping table missing columns: {', '.join(sorted(missing))}")
                continue
            row = dict(zip(header, fields))
            table_map[row["query"]] = {
                "renamed_to": row["renamed_to"],
                "needs_rc": row["needs_reverse_complement"].strip().lower() == "yes",
            }
    if not table_map:
        raise ValueError(f"Mapping table {mapping_table_path} is empty.")
    print(f"  Loaded {len(table_map)} entries from mapping table")

    # Build a HAP-agnostic lookup: strip _HAPJ from table keys so we can
    # match SUPER_1_HAP2 against a table that was built from SUPER_1_HAP1.
    base_table_map: Dict[str, Dict] = {}
    for q, info in table_map.items():
        if is_unloc_contig(q):
            continue
        base_q = strip_hap_suffix(q)
        base_table_map[base_q] = {
            "base_renamed_to": strip_hap_suffix(info["renamed_to"]),
            "needs_rc": info["needs_rc"],
        }

    # Same for unloc parent resolution
    parent_new_suffix = {
        strip_hap_suffix(q): _suffix_from_renamed(info["renamed_to"])
        for q, info in table_map.items()
        if not is_unloc_contig(q)
    }

    assignments: List[FinalChromosomeAssignment] = []
    for orig in (
        n for n in sequences
        if n.startswith(query_chromosome_prefix) and not is_unloc_contig(n)
    ):
        hap_match = _HAP_SUFFIX_RE.search(orig)
        hap_str = hap_match.group(0) if hap_match else ""
        base_orig = strip_hap_suffix(orig)

        if base_orig in base_table_map:
            entry = base_table_map[base_orig]
            suffix = _suffix_from_renamed(entry["base_renamed_to"])
            needs_rc = entry["needs_rc"]
        else:
            suffix = strip_hap_suffix(extract_chromosome_suffix(orig, query_chromosome_prefix))
            needs_rc = False
            print(f"  Warning: {orig} not in mapping table -- keeping original suffix")

        assignments.append(FinalChromosomeAssignment(
            original_name=orig,
            new_name=f"{output_prefix}{suffix}{hap_str}",
            new_suffix=suffix,
            needs_reverse_complement=needs_rc,
            is_sex_chromosome=is_sex_chromosome_suffix(suffix),
        ))

    unloc_mappings: List[UnlocMapping] = []
    for contig in (
        n for n in sequences
        if n.startswith(query_chromosome_prefix) and is_unloc_contig(n)
    ):
        parent, unloc_num = parse_unloc_name(contig)
        base_parent = strip_hap_suffix(parent)
        if base_parent not in parent_new_suffix:
            print(f"  Warning: parent '{parent}' for '{contig}' not in mapping table -- skipping")
            continue
        unloc_mappings.append(UnlocMapping(
            contig_name=contig, parent_chromosome=parent,
            unloc_number=unloc_num, needs_reverse_complement=False,
        ))

    print(f"  Built {len(assignments)} assignments and {len(unloc_mappings)} unloc mappings")
    return assignments, unloc_mappings
