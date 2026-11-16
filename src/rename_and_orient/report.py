"""Reports: text tables (returned as strings), mapping TSV and chromosome list files."""
import logging
from pathlib import Path

from .models import ChromosomeMapping, FinalChromosomeAssignment, UnlocMapping
from .names import _unloc_output_name, group_unlocs_by_parent, natural_sort_key

logger = logging.getLogger(__name__)


def format_assignment_table(assignments: list[FinalChromosomeAssignment]) -> str:
    """Format the final chromosome assignments as a text table.

    Args:
        assignments: Final chromosome assignments in output order.

    Returns:
        The table text, starting with a blank line and without a trailing newline.
    """
    lines = [
        "\nFinal chromosome assignments:",
        "-" * 70,
        f"{'Original':<20} {'New Name':<20} {'Suffix':<10} {'RC?':>5} {'Sex?':>5}",
        "-" * 70,
    ]
    for a in assignments:
        lines.append(
            f"{a.original_name:<20} {a.new_name:<20} {a.new_suffix:<10} "
            f"{'Yes' if a.needs_reverse_complement else 'No':>5} "
            f"{'Yes' if a.is_sex_chromosome else 'No':>5}"
        )
    return "\n".join(lines)


def format_mapping_summary(mappings: list[ChromosomeMapping]) -> str:
    """Format the chromosome mapping summary as a text table.

    Args:
        mappings: List of ChromosomeMapping objects.

    Returns:
        The table text, starting with a blank line and without a trailing newline.
    """
    lines = [
        "\nChromosome mapping summary:",
        "-" * 80,
        f"{'Query':<15} {'Target':<10} {'Coverage':>10} {'Plus':>12} {'Minus':>12} {'RC?':>5}",
        "-" * 80,
    ]
    for m in sorted(mappings, key=lambda x: natural_sort_key(x.target_name, "chr_")):
        lines.append(
            f"{m.query_name:<15} {m.target_name:<10} {m.coverage:>9.1%} "
            f"{m.plus_strand_length:>12,} {m.minus_strand_length:>12,} "
            f"{'Yes' if m.needs_reverse_complement else 'No':>5}"
        )
    return "\n".join(lines)


def save_mapping_tsv(
    mappings: list[ChromosomeMapping],
    assignments: list[FinalChromosomeAssignment],
    output_path: Path
) -> None:
    """Save chromosome mapping summary to TSV file.

    Args:
        mappings: List of ChromosomeMapping objects
        assignments: List of FinalChromosomeAssignment for renamed_to column
        output_path: Path to output TSV file
    """
    rename_lookup = {a.original_name: a.new_name for a in assignments}

    with open(output_path, 'w') as f:
        f.write("query\ttarget\trenamed_to\tquery_length\talignment_length\tcoverage\t"
                "plus_strand\tminus_strand\tneeds_reverse_complement\n")

        for m in sorted(mappings, key=lambda x: natural_sort_key(x.target_name, "chr_")):
            renamed_to = rename_lookup.get(m.query_name, m.query_name)
            f.write(f"{m.query_name}\t{m.target_name}\t{renamed_to}\t{m.query_length}\t"
                    f"{m.total_alignment_length}\t{m.coverage:.4f}\t"
                    f"{m.plus_strand_length}\t{m.minus_strand_length}\t"
                    f"{'yes' if m.needs_reverse_complement else 'no'}\n")

    logger.info("Mapping summary saved to: %s", output_path)


def write_chromosome_list(
    assignments: list[FinalChromosomeAssignment],
    unloc_mappings: list[UnlocMapping],
    output_path: Path,
    output_prefix: str = "SUPER_"
) -> None:
    """Write chromosome list CSV file.

    Format: name,suffix,yes/no (yes for main chr, no for unloc)

    Args:
        assignments: List of FinalChromosomeAssignment (sorted for output)
        unloc_mappings: List of UnlocMapping for unloc contigs
        output_path: Path to output CSV file
        output_prefix: Prefix for output chromosome names
    """
    parent_suffix_lookup = {a.original_name: a.new_suffix for a in assignments}

    unloc_by_parent = group_unlocs_by_parent(unloc_mappings)

    with open(output_path, 'w') as f:
        for a in assignments:
            f.write(f"{a.new_name},{a.new_suffix},yes\n")

            for unloc in unloc_by_parent.get(a.original_name, []):
                new_name, new_parent_suffix = _unloc_output_name(
                    unloc, parent_suffix_lookup, output_prefix
                )
                f.write(f"{new_name},{new_parent_suffix},no\n")

    logger.info("Chromosome list written to: %s", output_path)

