"""Load pre-built mapping table for second-haplotype renaming."""
import logging
from collections.abc import Collection
from pathlib import Path
from typing import Any

from .models import FinalChromosomeAssignment, UnlocMapping
from .names import (
    _HAP_SUFFIX_RE,
    extract_chromosome_suffix,
    is_sex_chromosome_suffix,
    is_unloc_contig,
    parse_unloc_name,
    strip_hap_suffix,
)

logger = logging.getLogger(__name__)


def load_mapping_table_assignments(
    mapping_table_path: Path,
    sequences: Collection[str],
    query_chromosome_prefix: str = "SUPER_",
    output_prefix: str = "SUPER_",
) -> tuple[list[FinalChromosomeAssignment], list[UnlocMapping]]:
    """Build assignments and unloc mappings from a pre-built mapping TSV.

    The TSV is produced by a previous run on haplotype 1. Columns used: query,
    renamed_to, needs_reverse_complement. Unlocs present only in this
    haplotype are mapped via their parent chromosome; they are never RC'd.
    """
    table_map: dict[str, dict[str, Any]] = {}
    with open(mapping_table_path) as fh:
        header = None
        for lineno, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if header is None:
                header = fields
                missing = {"query", "renamed_to", "needs_reverse_complement"} - set(header)
                if missing:
                    raise ValueError(f"Mapping table missing columns: {', '.join(sorted(missing))}")
                continue
            if len(fields) < len(header):
                raise ValueError(
                    f"Mapping table {mapping_table_path}, line {lineno}: expected "
                    f"{len(header)} tab-separated fields, got {len(fields)} "
                    f"(spaces instead of tabs?)"
                )
            row = dict(zip(header, fields))
            table_map[row["query"]] = {
                "renamed_to": row["renamed_to"],
                "needs_rc": row["needs_reverse_complement"].strip().lower() == "yes",
            }
    if not table_map:
        raise ValueError(f"Mapping table {mapping_table_path} is empty.")
    logger.info("Loaded %d entries from mapping table", len(table_map))

    def _suffix_from_renamed(renamed_to: str) -> str:
        """Extract the bare chromosomal suffix (e.g. '1', 'X', 'MT') from renamed_to.

        Strips output_prefix (or query_chromosome_prefix, which the table's earlier run
        may have used) and any _HAPn tag.
        """
        base = strip_hap_suffix(renamed_to)
        for prefix in (output_prefix, query_chromosome_prefix):
            if base.startswith(prefix):
                return base[len(prefix):]
        # shortcut: unknown prefix is cut at the last '_', so it breaks if the suffix itself has '_'
        return base.rsplit("_", 1)[-1]

    parent_new_suffix = {
        q: _suffix_from_renamed(info["renamed_to"])
        for q, info in table_map.items() if not is_unloc_contig(q)
    }

    assignments: list[FinalChromosomeAssignment] = []
    for orig in (
        n for n in sequences
        if n.startswith(query_chromosome_prefix) and not is_unloc_contig(n)
    ):
        if orig not in table_map:
            suffix = strip_hap_suffix(extract_chromosome_suffix(orig, query_chromosome_prefix))
            logger.warning("%s not in mapping table -- keeping original suffix", orig)
        else:
            suffix = _suffix_from_renamed(table_map[orig]["renamed_to"])
        needs_rc = table_map[orig]["needs_rc"] if orig in table_map else False
        # The _HAPn tag belongs to the scaffold, not to the chromosome number:
        # keep it in the output name but never in new_suffix (which is sorted as int).
        hap_match = _HAP_SUFFIX_RE.search(orig)
        hap_str = hap_match.group(0) if hap_match else ""
        assignments.append(FinalChromosomeAssignment(
            original_name=orig,
            new_name=f"{output_prefix}{suffix}{hap_str}",
            new_suffix=suffix,
            needs_reverse_complement=needs_rc,
            is_sex_chromosome=is_sex_chromosome_suffix(suffix),
        ))

    unloc_mappings: list[UnlocMapping] = []
    for contig in (
        n for n in sequences
        if n.startswith(query_chromosome_prefix) and is_unloc_contig(n)
    ):
        parent, unloc_num = parse_unloc_name(contig)
        if parent not in parent_new_suffix:
            logger.warning("parent '%s' for '%s' not in mapping table -- skipping", parent, contig)
            continue
        unloc_mappings.append(UnlocMapping(
            contig_name=contig, parent_chromosome=parent,
            unloc_number=unloc_num, needs_reverse_complement=False,
        ))

    logger.info(
        "Built %d assignments and %d unloc mappings", len(assignments), len(unloc_mappings)
    )
    return assignments, unloc_mappings
