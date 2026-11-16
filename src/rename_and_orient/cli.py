"""Command-line interface entry point."""
__version__ = "1.2.4"

import argparse
import logging
import sys
from pathlib import Path

from .alignment import build_chromosome_mappings
from .assignments import (
    build_unloc_mappings,
    resolve_chromosome_assignments,
    sort_assignments_for_output,
)
from .fasta import FastaIndex, calculate_genome_length, index_fasta, write_fasta
from .mapping_table import load_mapping_table_assignments
from .models import FinalChromosomeAssignment, UnlocMapping
from .paf import filter_paf_records, parse_paf, validate_paf_fasta_consistency
from .plotting import MissingPlotDependency, plot_alignments_overview
from .report import (
    format_assignment_table,
    format_mapping_summary,
    save_mapping_tsv,
    write_chromosome_list,
)

# Hard-coded: __name__ is "__main__" under `python -m rename_and_orient.cli`.
PACKAGE_LOGGER_NAME = "rename_and_orient"
logger = logging.getLogger(f"{PACKAGE_LOGGER_NAME}.cli")


class GenomeLengthMismatchError(RuntimeError):
    """Raised when the output genome length differs from the input genome length."""


class _ConsoleFormatter(logging.Formatter):
    """Plain message for INFO and below; ``LEVEL: message`` for WARNING and above."""

    def __init__(self) -> None:
        super().__init__()
        self._plain = logging.Formatter("%(message)s")
        self._leveled = logging.Formatter("%(levelname)s: %(message)s")

    def format(self, record: logging.LogRecord) -> str:
        """Format a record according to its level."""
        fmt = self._leveled if record.levelno >= logging.WARNING else self._plain
        return fmt.format(record)


def _setup_logging(
    verbose: bool = False, quiet: bool = False, log_file: Path | None = None
) -> None:
    """Configure the package logger (not the root logger).

    Console output goes to stderr. Handlers added by a previous call are removed first.

    Args:
        verbose: Show DEBUG messages on the console.
        quiet: Show only WARNING and above on the console.
        log_file: Optional file that receives all messages (DEBUG and above) with timestamps.
    """
    pkg_logger = logging.getLogger(PACKAGE_LOGGER_NAME)
    for handler in list(pkg_logger.handlers):
        if getattr(handler, "_rao_cli", False):
            pkg_logger.removeHandler(handler)
            handler.close()

    pkg_logger.setLevel(logging.DEBUG)

    console = logging.StreamHandler(sys.stderr)
    console.setLevel(logging.DEBUG if verbose else logging.WARNING if quiet else logging.INFO)
    console.setFormatter(_ConsoleFormatter())
    console._rao_cli = True  # type: ignore[attr-defined]
    pkg_logger.addHandler(console)

    if log_file is not None:
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
        )
        file_handler._rao_cli = True  # type: ignore[attr-defined]
        pkg_logger.addHandler(file_handler)


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
        help="Prefix for output chromosome names "
             "(e.g., SUPER_, chr_, chr, or empty string for no prefix)"
    )

    parser.add_argument(
        "--reference-chromosome-prefix", "-r",
        type=str,
        default=None,
        help="Prefix for reference chromosome names in PAF target "
             "(e.g., chr_, chr, SUPER_, scaffold_). "
             "Auto-detected from PAF if not specified."
    )

    parser.add_argument(
        "--plot-alignments", "-P",
        action="store_true",
        default=False,
        help="Write <output-prefix>.alignments_overview.pdf to the output directory: "
             "a grid of every chromosome's PAF alignment blocks, page 1 as aligned, "
             "page 2 after the orientation decision. Requires matplotlib."
    )

    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Show debug messages on stderr"
    )
    verbosity.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Show only warnings and errors on stderr (result tables are still printed)"
    )

    parser.add_argument(
        "--log-file",
        type=Path,
        default=None,
        help="Also write all log messages (debug level, with timestamps) to this file"
    )

    args = parser.parse_args()

    if not args.fasta.exists():
        parser.error(f"FASTA file not found: {args.fasta}")
    if args.mapping_table is not None and not args.mapping_table.exists():
        parser.error(f"Mapping table not found: {args.mapping_table}")

    if args.output_prefix is None:
        args.output_prefix = args.fasta.stem

    return args


def _write_and_validate(
    args: argparse.Namespace,
    sequences: FastaIndex,
    sorted_assignments: list[FinalChromosomeAssignment],
    unloc_mappings: list[UnlocMapping],
    output_chromosome_prefix: str,
    mapping_tsv_path: Path | None = None,
) -> None:
    """Print assignment table (stdout), write FASTA + CSV, validate genome length.

    Args:
        args: Parsed command-line arguments.
        sequences: Index of the input FASTA (sequences are read lazily, one at a time).
        sorted_assignments: Final chromosome assignments in output order.
        unloc_mappings: Unlocalised contig mappings.
        output_chromosome_prefix: Prefix for output chromosome names.
        mapping_tsv_path: Path of the mapping summary TSV, if one was written.

    Raises:
        GenomeLengthMismatchError: If the output genome length differs from the input.
    """
    print(format_assignment_table(sorted_assignments))

    fasta_out = args.output_dir / f"{args.output_prefix}.fa"
    csv_out = args.output_dir / f"{args.output_prefix}.chromosome.list.csv"

    logger.info("Writing output files...")
    write_fasta(sequences, sorted_assignments, unloc_mappings, fasta_out, output_chromosome_prefix)
    write_chromosome_list(sorted_assignments, unloc_mappings, csv_out, output_chromosome_prefix)

    logger.info("Validating genome length...")
    in_len = sequences.total_length
    out_len = calculate_genome_length(fasta_out)
    if in_len == out_len:
        logger.info("OK: Genome length matches (%s bp)", f"{in_len:,}")
    else:
        raise GenomeLengthMismatchError(
            f"Genome length mismatch! Input: {in_len:,} bp, Output: {out_len:,} bp"
        )

    logger.info("Done! Output files:")
    logger.info("  - FASTA: %s", fasta_out)
    logger.info("  - Chromosome list: %s", csv_out)
    if mapping_tsv_path:
        logger.info("  - Mapping summary: %s", mapping_tsv_path)


def main() -> None:
    """Main entry point."""
    args = parse_args()
    _setup_logging(args.verbose, args.quiet, args.log_file)
    try:
        _run(args)
    except GenomeLengthMismatchError as exc:
        logger.error("%s", exc)
        sys.exit(1)


def _run(args: argparse.Namespace) -> None:
    """Run the renaming pipeline for parsed arguments.

    Args:
        args: Parsed command-line arguments.
    """
    logger.info("rename-and-orient version %s", __version__)

    qpfx = args.query_chromosome_prefix
    opfx = args.output_chromosome_prefix

    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Reading FASTA file: %s", args.fasta)
    with index_fasta(args.fasta) as sequences:
        logger.info("Found %d sequences", len(sequences))
        _run_with_index(args, sequences, qpfx, opfx)


def _run_with_index(
    args: argparse.Namespace, sequences: FastaIndex, qpfx: str, opfx: str
) -> None:
    """Run the pipeline given the FASTA index.

    Args:
        args: Parsed command-line arguments.
        sequences: Index of the input FASTA.
        qpfx: Input chromosome prefix.
        opfx: Output chromosome prefix.
    """
    # --- Mode A: mapping table (second haplotype, no alignment needed) ---
    if args.mapping_table is not None:
        if args.min_coverage != 0.5:
            logger.warning("--min-coverage is ignored in --mapping-table mode")
        if args.plot_alignments:
            logger.warning("--plot-alignments is ignored in --mapping-table mode")
        if args.reference_chromosome_prefix is not None:
            logger.warning("--reference-chromosome-prefix is ignored in --mapping-table mode")
        logger.info("Using pre-built mapping table: %s", args.mapping_table)
        logger.info("Input prefix: '%s' -> Output prefix: '%s'", qpfx, opfx)
        assignments, unloc_mappings = load_mapping_table_assignments(
            args.mapping_table, sequences, qpfx, opfx
        )
        _write_and_validate(args, sequences, sort_assignments_for_output(assignments),
                            unloc_mappings, opfx)
        return

    # --- Mode B: full PAF-based alignment ---
    logger.info("Parsing PAF file: %s", args.paf)
    paf_records = parse_paf(args.paf)
    logger.info("Found %d alignment records", len(paf_records))

    ref_prefix_arg = args.reference_chromosome_prefix
    filtered_records, ref_prefix = filter_paf_records(paf_records, qpfx, ref_prefix_arg)
    label = "provided" if ref_prefix_arg else "auto-detected"
    logger.info("Reference prefix (%s): '%s'", label, ref_prefix)
    logger.info(
        "After filtering (%s* -> %s*): %d records", qpfx, ref_prefix, len(filtered_records)
    )

    logger.info("Validating PAF/FASTA consistency...")
    is_valid, in_paf_not_fasta, in_fasta_not_paf = validate_paf_fasta_consistency(
        paf_records, sequences, qpfx
    )
    if in_paf_not_fasta:
        logger.warning("In PAF but not FASTA: %s", ", ".join(in_paf_not_fasta))
    if in_fasta_not_paf:
        logger.warning(
            "In FASTA but not PAF: %s (will keep original suffix)", ", ".join(in_fasta_not_paf)
        )
    if is_valid:
        logger.info("OK: All %s* chromosomes match", qpfx)

    logger.info("Building chromosome mappings (min coverage: %.0f%%)...", args.min_coverage * 100)
    mappings = build_chromosome_mappings(filtered_records, args.min_coverage, ref_prefix)
    logger.info("Successfully mapped %d chromosomes", len(mappings))

    if args.plot_alignments:
        logger.info("Generating alignment overview PDF...")
        try:
            plot_alignments_overview(
                filtered_records, mappings,
                args.output_dir / f"{args.output_prefix}.alignments_overview.pdf",
            )
        except MissingPlotDependency as exc:
            logger.warning("%s", exc)

    print(format_mapping_summary(mappings))

    logger.info("Resolving chromosome assignments...")
    logger.info("Input prefix: '%s' -> Output prefix: '%s'", qpfx, opfx)
    assignments, _ = resolve_chromosome_assignments(mappings, sequences, qpfx, opfx)

    unloc_mappings = build_unloc_mappings(sequences, mappings)
    if unloc_mappings:
        logger.info("Found %d unlocalized contigs", len(unloc_mappings))
    for unloc in unloc_mappings:
        unloc.needs_reverse_complement = False

    sorted_assignments = sort_assignments_for_output(assignments)
    mapping_tsv_path = args.output_dir / f"{args.output_prefix}.mapping.tsv"
    save_mapping_tsv(mappings, sorted_assignments, mapping_tsv_path)

    _write_and_validate(args, sequences, sorted_assignments, unloc_mappings, opfx, mapping_tsv_path)


if __name__ == "__main__":
    main()
