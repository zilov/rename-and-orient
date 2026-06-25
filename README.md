# rename-and-orient

Rename and orient chromosomes in a FASTA file based on alignment to a reference genome (PAF format).

Developed for genome curation workflows in the [Darwin Tree of Life](https://www.darwintreeoflife.org/) project (GRIT team, Wellcome Sanger Institute).

## Features

- Renames scaffolds (`SUPER_N`) to chromosome names based on best PAF alignment to reference
- Reverse-complements chromosomes to match reference orientation (Pearson correlation method)
- Handles sex chromosomes (W, Z, X, Y, Z1/Z2, etc.), unlocalized contigs (`_unloc_`), and HAP-suffixed scaffolds
- Resolves assignment conflicts (multiple query chromosomes mapping to the same target)
- `--mapping-table` mode: rename a second haplotype using the mapping produced from the first (no re-alignment needed)
- Optional alignment scatter plots (`--plot-alignments`, requires matplotlib)

## Installation

```bash
pip install rename-and-orient
# or with uv:
uv tool install rename-and-orient
```

## Usage

### PAF mode (standard)

```bash
rename-and-orient \
  --fasta hap1.primary.curated.fa \
  --paf hap1_vs_reference.paf \
  --output-dir results/ \
  --output-prefix ilMysSpe1.hap1.1
```

### Mapping-table mode (second haplotype)

```bash
rename-and-orient \
  --fasta hap2.primary.curated.fa \
  --mapping-table results/ilMysSpe1.hap1.1.mapping.tsv \
  --output-dir results_hap2/ \
  --output-prefix ilMysSpe1.hap2.1
```

## Options

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--fasta` | `-f` | required | Input FASTA (gzip supported) |
| `--paf` | `-p` | required* | PAF alignment file |
| `--mapping-table` | `-mt` | required* | Pre-built mapping TSV (mutually exclusive with `--paf`) |
| `--output-dir` | `-d` | `./rename_and_orient` | Output directory |
| `--output-prefix` | `-o` | FASTA stem | Output file prefix |
| `--min-coverage` | `-c` | `0.5` | Minimum alignment coverage to rename |
| `--query-chromosome-prefix` | `-q` | `SUPER_` | Prefix for query scaffolds |
| `--output-chromosome-prefix` | `-x` | `SUPER_` | Prefix for output chromosome names |
| `--reference-chromosome-prefix` | `-r` | auto | Reference chromosome prefix in PAF |
| `--plot-alignments` | `-P` | off | Generate alignment scatter plots (requires matplotlib) |

## Output files

| File | Description |
|------|-------------|
| `<prefix>.fa` | Renamed and oriented FASTA |
| `<prefix>.chromosome.list.csv` | Chromosome list CSV (`name,suffix,yes/no`) |
| `<prefix>.mapping.tsv` | Mapping summary TSV (PAF mode only) |
| `plots/` | Alignment scatter plots (if `--plot-alignments`) |

## Development

```bash
git clone https://github.com/zilov/rename-and-orient
cd rename-and-orient
uv sync
uv run pytest tests/
```
