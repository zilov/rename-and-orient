"""Load pre-built mapping table for second-haplotype renaming."""
import re
from pathlib import Path
from typing import Dict, List, Tuple

from .models import FinalChromosomeAssignment, UnlocMapping
from .names import extract_chromosome_suffix, is_sex_chromosome_suffix, is_unloc_contig, parse_unloc_name


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

    def _suffix_from_renamed(renamed_to: str) -> str:
        """Extract the bare chromosomal suffix (e.g. '1', 'X') from renamed_to
        regardless of which prefix was used in the original run."""
        m = re.search(r'([A-Z]\d*|\d+)$', renamed_to, re.IGNORECASE)
        return m.group(1) if m else renamed_to

    parent_new_suffix = {
        q: _suffix_from_renamed(info["renamed_to"])
        for q, info in table_map.items() if not is_unloc_contig(q)
    }

    assignments: List[FinalChromosomeAssignment] = []
    for orig in (
        n for n in sequences
        if n.startswith(query_chromosome_prefix) and not is_unloc_contig(n)
    ):
        if orig not in table_map:
            suffix = extract_chromosome_suffix(orig, query_chromosome_prefix)
            print(f"  Warning: {orig} not in mapping table -- keeping original suffix")
        else:
            suffix = _suffix_from_renamed(table_map[orig]["renamed_to"])
        needs_rc = table_map[orig]["needs_rc"] if orig in table_map else False
        assignments.append(FinalChromosomeAssignment(
            original_name=orig,
            new_name=f"{output_prefix}{suffix}",
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
        if parent not in parent_new_suffix:
            print(f"  Warning: parent '{parent}' for '{contig}' not in mapping table -- skipping")
            continue
        unloc_mappings.append(UnlocMapping(
            contig_name=contig, parent_chromosome=parent,
            unloc_number=unloc_num, needs_reverse_complement=False,
        ))

    print(f"  Built {len(assignments)} assignments and {len(unloc_mappings)} unloc mappings")
    return assignments, unloc_mappings
