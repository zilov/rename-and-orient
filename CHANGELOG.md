# Changelog

All notable changes to this project are documented in this file.

Format based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

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
