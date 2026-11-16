# Ponytail audit — 2026-11-16

**What this repo does:** `rename-and-orient` is a CLI tool for Darwin Tree of Life / GRIT. It reads a FASTA file and a PAF alignment to a reference, renames `SUPER_N` chromosomes to match the reference, and reverse-complements the ones that lie the wrong way round. The `--mapping-table` mode applies a ready-made table to a second haplotype.

**Assumed load:** one person runs the script on one genome (hundreds of MB to a few GB). No parallel processes, no state shared between processes. Tests: 147 passed, 9 skipped.

## Must fix

### 1. Running the tool can destroy the input FASTA (`cli.py:L227-228`, `fasta.py:L274`)
- **What this is:** The output file is named `<output-prefix>.fa`. The default prefix is the input file name without its extension.
- **Problem:** For `genome.fa` with `-d .`, the output path equals the input path. `open(output, 'w')` truncates the file while `FastaIndex` is still reading from it. Verified: the 18-byte input became headers only, the program exited with "Genome length mismatch", and the original data was already gone.
- **Fix:** In `_run`, before writing, compare `args.output_dir.resolve() / f"{prefix}.fa"` with `args.fasta.resolve()`. If they are equal, stop with a clear error. About 3 lines.
- **If we skip it:** One day a curator loses the only copy of an assembly.

### 2. Crash on an "unusual" chromosome name (`assignments.py:L225`, `mapping_table.py:L80-86`)
- **What this is:** Before writing files, chromosomes are sorted by number with `int(x.new_suffix)`.
- **Problem:** An unmapped scaffold `SUPER_1A` raises `ValueError: '1A'` (verified). In mapping-table mode, the regex turns `SUPER_MT` into the suffix `T`, with the same crash. The program dies with a traceback after all the work and writes nothing.
- **Fix:** Treat a non-numeric suffix as "not an autosome": in `sort_assignments_for_output` use the key `(0, int(s)) if s.isdigit() else (1, s)`. In `_suffix_from_renamed`, take the suffix by `output_prefix` instead of guessing with a regex.
- **If we skip it:** Any genome with a non-standard scaffold name cannot be processed without hand-editing the input.

### 3. Release to PyPI without running tests (`.github/workflows/release.yml`)
- **What this is:** Every push to `main` with a new version is published to PyPI.
- **Problem:** The workflow has no `pytest` step. Also, `matplotlib` is not in the dev group, so the 9 plotting tests are skipped both on developer machines and in CI.
- **Fix:** Add `uv run --extra plots pytest` before `uv build`. Add `matplotlib` to `[dependency-groups].dev`.
- **If we skip it:** A broken version goes out to everyone who runs `pip install -U`.

## Should fix

### 4. A bad PAF gives a cryptic error, and short lines are skipped silently (`paf.py:L52-62`)
- **What this is:** The parser reads the 12 PAF columns.
- **Problem:** A non-numeric field gives a bare `ValueError: 'X'` with no line number. Lines with fewer than 12 columns are dropped without a message, so a truncated file silently gives an incomplete result.
- **Fix:** Wrap the line parsing in `try/except ValueError` and raise an error with the line number. For short lines, keep a counter and write one `logger.warning`.
- **If we skip it:** An incomplete PAF leads to wrong or missing renames with no warning.

### 5. Version is written in two places (`pyproject.toml:L3`, `cli.py:L2`)
- **What this is:** The release takes its version from `pyproject.toml`, while `--version` shows `cli.__version__`.
- **Problem:** Both are edited by hand. If one is forgotten, release `1.3.0` will print `--version 1.2.4`.
- **Fix:** `__version__ = importlib.metadata.version("rename-and-orient")` and delete the constant.
- **If we skip it:** Logs and `--version` show the wrong version, and it is unclear what was actually run.

### 6. No range check for `--min-coverage` (`cli.py:L132-137`)
- **What this is:** The coverage threshold accepts any number.
- **Problem:** `-c 50` (meaning 50%) makes the threshold impossible: no chromosome is renamed, yet the program finishes "successfully".
- **Fix:** `if not 0 <= args.min_coverage <= 1: parser.error(...)`.
- **If we skip it:** A typo gives a "successful" run with an empty result.

### 7. Only one error type is caught (`cli.py:L262-267`)
- **What this is:** `main` catches only `GenomeLengthMismatchError`.
- **Problem:** Other `ValueError`s (empty table, duplicate FASTA name, no prefix in PAF) escape as a traceback, although their messages are already written for humans.
- **Fix:** Catch `(GenomeLengthMismatchError, ValueError)`, call `logger.error` and `sys.exit(1)`.
- **If we skip it:** The user sees a page of traceback instead of one clear line.

### 8. Dead and duplicated code (about -70 lines)
- **What this is:** Leftovers from the refactoring.
- **Problem:**
  - `io.py` (32 lines) is a deprecated shim; nothing in the repo imports it except the test of the shim itself.
  - The re-exports at the end of `alignment.py` (L186-192) exist only for old imports.
  - `report.print_mapping_summary` is a one-call wrapper.
  - `PAFRecord.alignment_block_length` is never used.
  - The second return value of `resolve_chromosome_assignments` (`rc_lookup`) is thrown away by the CLI.
  - Grouping unlocs by parent is written twice (`fasta.py:L249-254`, `report.py:L109-113`).
- **Fix:** Delete `io.py` and its test, the `alignment.py` re-exports, `print_mapping_summary` and `alignment_block_length`. Return only the list from `resolve_chromosome_assignments`. Move the grouping into `names.group_unlocs_by_parent`. If outside code depends on the shim, keep `io.py` for one more version.
- **If we skip it:** Every one of these pieces must still be read, tested and maintained for nothing.

### 9. Two plotting paths (`plotting.py`, about 140 lines)
- **What this is:** The CLI uses only `plot_alignments_overview`. A separate per-chromosome branch (`render_alignment_figure` / `save_figure` / `plot_chromosome_alignments`, one file per chromosome) is kept "as a library function".
- **Problem:** It is an interface with no call from real code (tests only). It doubles the amount of plotting code to maintain.
- **Fix:** If there are no outside users, delete the per-chromosome branch and its tests. This is the author's decision.
- **If we skip it:** Two renderers are maintained where one is enough.

## Nice to have

### 10. Small things
- `"_unloc_" in name` is written directly in 5 places next to the ready-made `is_unloc_contig` (`assignments.py`, `paf.py`).
- `merge_intervals(strands['+'\n ])` in `alignment.py:L77-79` is split over two lines by mistake.
- About 3 smaller findings were left out.

## Verdict
Fix 1 first (protect the input file from being overwritten), then 2 and 3. The core logic (Pearson orientation, interval merging, lazy FASTA index) looks clean and is covered by tests. I found no errors in the target-choice or orientation logic itself.

**Lean:** -70 lines, -0 dependencies possible (there are no dependencies; `matplotlib` is optional). Removing the plotting branch from finding 9 would save about -140 more lines.

## Not checked
- `plotting.py` from L120 (panel layout, colours).
- Tests were read only by name and by run result.
- Real large genomes and gzip input of hundreds of GB were not run.
- `TODO/tol_audit_2026-11-16.md` was not read, to avoid relying on someone else's conclusions.
