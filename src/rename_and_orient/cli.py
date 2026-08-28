"""Command-line interface entry point."""
__version__ = "1.2.2"

import argparse
import sys
from pathlib import Path

from .alignment import (
    build_chromosome_mappings,
    filter_paf_records,
    plot_chromosome_alignments,
    validate_paf_fasta_consistency,
)
from .assignments import build_unloc_mappings, resolve_chromosome_assignments, sort_assignments_for_output
from .io import (
    calculate_genome_length,
    print_mapping_summary,
    read_fasta,
    save_mapping_tsv,
    write_chromosome_list,
    write_fasta,
)
from .mapping_table import load_mapping_table_assignments


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Rename and orient chromosomes based on reference alignment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--version", "-v",
        action="version",
        version=f"%(prog)s {__version__}"
    )

    parser.add_argument(
        "--fasta", "-f",
        required=True,
        type=Path,
        help="Input FASTA file (can be gzipped)"
    )

    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "--paf", "-p",
        type=Path,
        default=None,
        help="PAF file with alignment to reference"
    )
    source_group.add_argument(
        "--mapping-table", "-mt",
        type=Path,
        default=None,
        dest="mapping_table",
        help="Pre-built mapping TSV (output of a previous run) to rename/orient a second "
             "haplotype without re-running alignment. Columns used: query, renamed_to, "
             "needs_reverse_complement."
    )

    parser.add_argument(
        "--output-dir", "-d",
        type=Path,
        default=Path("./rename_and_orient"),
        help="Output directory for generated files (will be created if it doesn't exist)"
    )

    parser.add_argument(
        "--output-prefix", "-o",
        type=str,
        help="Prefix for output file names (default: derived from input FASTA file name)"
    )

    parser.add_argument(
        "--min-coverage", "-c",
        type=float,
        default=0.5,
        help="Minimum coverage threshold for renaming (0.0-1.0)"
    )

    parser.add_argument(
        "--query-chromosome-prefix", "-q",
        type=str,
        default="SUPER_",
        help="Prefix for query chromosome names in input FASTA (e.g., SUPER_, scaffold_, contig_)"
    )

    parser.add_argument(
        "--output-chromosome-prefix", "-x",
        type=str,
        default="SUPER_",
        help="Prefix for output chromosome names (e.g., SUPER_, chr_, chr, or empty string for no prefix)"
    )

    parser.add_argument(
        "--reference-chromosome-prefix", "-r",
        type=str,
        default=None,
        help="Prefix for reference chromosome names in PAF target (e.g., chr_, chr, SUPER_, scaffold_). "
             "Auto-detected from PAF if not specified."
    )

    parser.add_argument(
        "--plot-alignments", "-P",
        action="store_true",
        default=False,
        help="Generate scatter plots of PAF alignment blocks for each chromosome "
             "(saved to <output-dir>/plots/). Requires matplotlib."
    )

    args = parser.parse_args()

    if args.output_prefix is None:
        args.output_prefix = args.fasta.stem

    return args


def _write_and_validate(
    args: argparse.Namespace,
    sequences: dict,
    sorted_assignments: list,
    unloc_mappings: list,
    output_chromosome_prefix: str,
    mapping_tsv_path: Path = None,
) -> None:
    """Print assignment table, write FASTA + CSV, validate genome length."""
    print("\nFinal chromosome assignments:")
    print("-" * 70)
    print(f"{'Original':<20} {'New Name':<20} {'Suffix':<10} {'RC?':>5} {'Sex?':>5}")
    print("-" * 70)
    for a in sorted_assignments:
        print(f"{a.original_name:<20} {a.new_name:<20} {a.new_suffix:<10} "
              f"{'Yes' if a.needs_reverse_complement else 'No':>5} "
              f"{'Yes' if a.is_sex_chromosome else 'No':>5}")

    fasta_out = args.output_dir / f"{args.output_prefix}.fa"
    csv_out = args.output_dir / f"{args.output_prefix}.chromosome.list.csv"

    print("\nWriting output files...")
    write_fasta(sequences, sorted_assignments, unloc_mappings, fasta_out, output_chromosome_prefix)
    write_chromosome_list(sorted_assignments, unloc_mappings, csv_out, output_chromosome_prefix)

    print("\nValidating genome length...")
    in_len = sum(len(s) for s in sequences.values())
    out_len = calculate_genome_length(fasta_out)
    if in_len == out_len:
        print(f"  OK: Genome length matches ({in_len:,} bp)")
    else:
        print(f"  ERROR: Genome length mismatch! Input: {in_len:,} bp, Output: {out_len:,} bp",
              file=sys.stderr)
        sys.exit(1)

    print(f"\nDone! Output files:")
    print(f"  - FASTA: {fasta_out}")
    print(f"  - Chromosome list: {csv_out}")
    if mapping_tsv_path:
        print(f"  - Mapping summary: {mapping_tsv_path}")


def main():
    """Main entry point."""
    args = parse_args()

    print(f"rename-and-orient version {__version__}")

    qpfx = args.query_chromosome_prefix
    opfx = args.output_chromosome_prefix

    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not args.fasta.exists():
        print(f"Error: FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading FASTA file: {args.fasta}")
    sequences = read_fasta(args.fasta)
    print(f"  Found {len(sequences)} sequences")

    # --- Mode A: mapping table (second haplotype, no alignment needed) ---
    if args.mapping_table is not None:
        if not args.mapping_table.exists():
            print(f"Error: Mapping table not found: {args.mapping_table}", file=sys.stderr)
            sys.exit(1)
        if args.min_coverage != 0.5:
            print("  Warning: --min-coverage is ignored in --mapping-table mode", file=sys.stderr)
        if args.plot_alignments:
            print("  Warning: --plot-alignments is ignored in --mapping-table mode", file=sys.stderr)
        if args.reference_chromosome_prefix is not None:
            print("  Warning: --reference-chromosome-prefix is ignored in --mapping-table mode",
                  file=sys.stderr)
        print(f"\nUsing pre-built mapping table: {args.mapping_table}")
        print(f"  Input prefix: '{qpfx}' -> Output prefix: '{opfx}'")
        assignments, unloc_mappings = load_mapping_table_assignments(
            args.mapping_table, sequences, qpfx, opfx
        )
        _write_and_validate(args, sequences, sort_assignments_for_output(assignments),
                            unloc_mappings, opfx)
        return

    # --- Mode B: full PAF-based alignment ---
    print(f"Parsing PAF file: {args.paf}")
    paf_records = read_paf_for_cli(args.paf)
    print(f"  Found {len(paf_records)} alignment records")

    ref_prefix_arg = args.reference_chromosome_prefix
    filtered_records, ref_prefix = filter_paf_records(paf_records, qpfx, ref_prefix_arg)
    label = "provided" if ref_prefix_arg else "auto-detected"
    print(f"  Reference prefix ({label}): '{ref_prefix}'")
    print(f"  After filtering ({qpfx}* -> {ref_prefix}*): {len(filtered_records)} records")

    print("\nValidating PAF/FASTA consistency...")
    is_valid, in_paf_not_fasta, in_fasta_not_paf = validate_paf_fasta_consistency(
        paf_records, sequences, qpfx
    )
    if in_paf_not_fasta:
        print(f"  WARNING: In PAF but not FASTA: {', '.join(in_paf_not_fasta)}")
    if in_fasta_not_paf:
        print(f"  WARNING: In FASTA but not PAF: {', '.join(in_fasta_not_paf)} (will keep original suffix)")
    if is_valid:
        print(f"  OK: All {qpfx}* chromosomes match")

    print(f"\nBuilding chromosome mappings (min coverage: {args.min_coverage:.0%})...")
    mappings = build_chromosome_mappings(filtered_records, args.min_coverage, ref_prefix)
    print(f"  Successfully mapped {len(mappings)} chromosomes")

    if args.plot_alignments:
        print("\nGenerating alignment scatter plots...")
        plot_chromosome_alignments(filtered_records, mappings, args.output_dir)

    print_mapping_summary(mappings)

    print("\nResolving chromosome assignments...")
    print(f"  Input prefix: '{qpfx}' -> Output prefix: '{opfx}'")
    assignments, _ = resolve_chromosome_assignments(mappings, sequences, qpfx, opfx)

    unloc_mappings = build_unloc_mappings(sequences, mappings)
    if unloc_mappings:
        print(f"  Found {len(unloc_mappings)} unlocalized contigs")
    for unloc in unloc_mappings:
        unloc.needs_reverse_complement = False

    sorted_assignments = sort_assignments_for_output(assignments)
    mapping_tsv_path = args.output_dir / f"{args.output_prefix}.mapping.tsv"
    save_mapping_tsv(mappings, sorted_assignments, mapping_tsv_path)

    _write_and_validate(args, sequences, sorted_assignments, unloc_mappings, opfx, mapping_tsv_path)


def read_paf_for_cli(paf_path: Path):
    """Read PAF file for CLI usage (thin wrapper around io.parse_paf)."""
    from .io import parse_paf
    return parse_paf(paf_path)


if __name__ == "__main__":
    main()
