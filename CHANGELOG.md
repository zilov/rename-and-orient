# Changelog

All notable changes to this project are documented in this file.

Format based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Changed

- Progress messages and warnings now go through the standard `logging` module instead of `print`. They are written to **stderr** (previously some went to stdout); the result tables (mapping summary and final chromosome assignments) stay on **stdout**. The library never configures logging; a `NullHandler` is installed on the `rename_and_orient` logger.
- Table formatting moved to the new `rename_and_orient.report` module (`format_assignment_table`, `format_mapping_summary`). `io.print_mapping_summary` is kept as a thin wrapper.
- A genome length mismatch is raised as `GenomeLengthMismatchError` and reported by `main()` (exit code 1, as before).
- A missing `--fasta` / `--mapping-table` file is now reported by the argument parser (exit code 2 instead of 1).
- `reverse_complement` uses `str.translate` and now complements all IUPAC ambiguity codes (R<->Y, K<->M, H<->D, B<->V, S, W, N; both cases) instead of passing them through unchanged. Output for ACGTN-only sequences is identical to before; sequences containing ambiguity codes now differ. Table taken from agp-tpf-utils (MIT).
- The input FASTA is no longer loaded into memory: `fasta.index_fasta()` indexes it in one pass and records are read and written one at a time, so peak memory is roughly that of the largest single record (2-3x its size) instead of the whole genome (twice, during the length check). Output files are byte-identical to before. `write_fasta` accepts either a dict or a `FastaIndex`; the name-only functions (`validate_paf_fasta_consistency`, `resolve_chromosome_assignments`, `build_unloc_mappings`, `load_mapping_table_assignments`) accept any collection of names. `calculate_genome_length` streams the file.
- Duplicate sequence names in the input FASTA now raise a `ValueError` (previously the last record silently won). This also applies to `read_fasta`.
- A gzipped input FASTA is decompressed once to a temporary file in `$TMPDIR` (removed when the run ends), because gzip does not support cheap random access. Plain FASTA files are read in place.

### Added
- `--plot-alignments` now writes a single `<output-prefix>.alignments_overview.pdf` next to the other outputs instead of per-chromosome images in `plots/`: a two-page PDF with a grid of every chromosome pair's alignment blocks, page 1 as aligned ("before"), page 2 after the orientation decision (target axis mirrored where RC is decided). Panels are ordered by reference chromosome and rasterized to keep the file small. `plotting.plot_chromosome_alignments` remains available as a library function.

- New modules `fasta`, `paf`, `orientation`, `plotting`; `names` gained `natural_sort_key`, `report` gained `save_mapping_tsv`, `write_chromosome_list` and `print_mapping_summary`. The public API in `rename_and_orient.__init__` is unchanged.
- `plotting.PlotData`, `prepare_plot_data`, `render_alignment_figure`, `save_figure` and `MissingPlotDependency` (an `ImportError`) for plotting with or without matplotlib.
- `--verbose` (debug messages), `--quiet` (warnings and errors only; result tables are still printed) and `--log-file PATH` (all messages at debug level, with timestamps). No short flags: `-v` is `--version` and `-q` is `--query-chromosome-prefix`.

- Module split: FASTA, PAF, orientation, reporting and plotting code moved out of `io.py` and `alignment.py` into dedicated modules. `alignment.py` still re-exports the moved names.
- `plot_chromosome_alignments` is now a thin wrapper around `prepare_plot_data` / `render_alignment_figure` / `save_figure`. It uses `matplotlib.figure.Figure` with the Agg canvas (no `pyplot`, no global `matplotlib.use`), accepts `fmt` and returns the list of written paths. A missing matplotlib raises `MissingPlotDependency`; the CLI logs a warning and continues without plots.

- `fasta.index_fasta`, `fasta.FastaIndex` (lazy read-only mapping name -> sequence with `names`, `length()`, `total_length`, `read_sequence()`; a context manager).

### Deprecated

- `rename_and_orient.io` is a deprecated shim that re-exports the old names and emits a `DeprecationWarning` on import. Import from `rename_and_orient.fasta`, `.paf`, `.report` or `.names` instead.

### Fixed

- Alignment plots use the real target length (`target_length` from the PAF) to mirror the right-hand panel, instead of `max(target_end)`.

## [1.2.4] - 2026-10-26

### Fixed

- `--mapping-table` ignores blank lines. A trailing empty line in the TSV used to crash with `KeyError: 'renamed_to'`.
- A mapping-table row with too few tab-separated fields (e.g. spaces instead of tabs) now raises a `ValueError` naming the file and line number, instead of a bare `KeyError`.
- Windows (`\r\n`) line endings in the mapping table are handled.

## [1.2.3] - 2026-09-22

### Fixed

- `--mapping-table` on a haplotype-tagged assembly no longer crashes in `sort_assignments_for_output` (`invalid literal for int() ... '1_HAP1'`). `new_suffix` always holds the bare chromosome suffix; the `_HAPn` tag is kept in `new_name`.
- `_suffix_from_renamed()` no longer returns `P1` for `SUPER_21_HAP1` / `SUPER_Z_HAP1`, which silently produced names like `SUPER_P1`.
- Unloc output names keep the `_HAPn` tag (`SUPER_21_unloc_1_HAP1`).

### Documentation

- README documents the minimal mapping-table format; only scaffolds to be renamed or flipped need to be listed.

## [1.2.2] - 2026-08-28

### Added

- The tool prints its version at startup.

## [1.2.1] - 2026-08-28

### Fixed

- `parse_unloc_name` handles haplotype-tagged unloc contigs such as `SUPER_1_unloc_2_HAP1`, which previously crashed the unloc-number parsing and mismatched the parent chromosome lookup.
