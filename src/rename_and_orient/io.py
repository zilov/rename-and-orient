"""I/O functions: reading FASTA/PAF, writing outputs."""
import gzip
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from .models import ChromosomeMapping, FinalChromosomeAssignment, PAFRecord, UnlocMapping
from .names import extract_chromosome_suffix, is_sex_chromosome_suffix, is_unloc_contig


def read_fasta(fasta_path: Path) -> Dict[str, str]:
    """
    Read FASTA file (supports gzip compression).

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Dictionary mapping sequence names to sequences
    """
    sequences = {}
    current_name = None
    current_seq = []

    open_func = gzip.open if str(fasta_path).endswith('.gz') else open
    mode = 'rt' if str(fasta_path).endswith('.gz') else 'r'

    with open_func(fasta_path, mode) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                if current_name is not None:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_name is not None:
            sequences[current_name] = ''.join(current_seq)

    return sequences


def parse_paf(paf_path: Path) -> List[PAFRecord]:
    """
    Parse PAF file.

    PAF format (tab-separated):
    1. Query sequence name
    2. Query sequence length
    3. Query start (0-based)
    4. Query end (0-based, open)
    5. Strand ('+' or '-')
    6. Target sequence name
    7. Target sequence length
    8. Target start (0-based)
    9. Target end (0-based, open)
    10. Number of matching bases
    11. Alignment block length
    12. Mapping quality (0-255)

    Args:
        paf_path: Path to PAF file

    Returns:
        List of PAFRecord objects
    """
    records = []

    with open(paf_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')
            if len(fields) < 12:
                continue

            record = PAFRecord(
                query_name=fields[0],
                query_length=int(fields[1]),
                query_start=int(fields[2]),
                query_end=int(fields[3]),
                strand=fields[4],
                target_name=fields[5],
                target_length=int(fields[6]),
                target_start=int(fields[7]),
                target_end=int(fields[8]),
                num_matches=int(fields[9]),
                alignment_length=int(fields[10]),
                mapping_quality=int(fields[11])
            )
            records.append(record)

    return records


def reverse_complement(seq: str) -> str:
    """
    Return reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        Reverse complement sequence
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}

    result = []
    for base in reversed(seq):
        result.append(complement.get(base, base))

    return ''.join(result)


def natural_sort_key(name: str, prefix: str = "SUPER_") -> Tuple:
    """
    Generate a sort key for natural sorting of chromosome names.

    Sorts autosomes numerically (1, 2, ... 10, 11, ...)
    and sex chromosomes alphabetically after autosomes.

    Args:
        name: Chromosome name (e.g., "SUPER_1", "SUPER_W")
        prefix: Prefix to strip (default "SUPER_")

    Returns:
        Tuple for sorting: (is_sex_chr, numeric_or_alpha_key)
    """
    suffix = extract_chromosome_suffix(name, prefix)

    if is_sex_chromosome_suffix(suffix):
        return (1, suffix)
    else:
        try:
            return (0, int(suffix))
        except ValueError:
            return (0, float('inf'))


def print_mapping_summary(mappings: List[ChromosomeMapping]) -> None:
    """
    Print summary of chromosome mappings.

    Args:
        mappings: List of ChromosomeMapping objects
    """
    print("\nChromosome mapping summary:")
    print("-" * 80)
    print(f"{'Query':<15} {'Target':<10} {'Coverage':>10} {'Plus':>12} {'Minus':>12} {'RC?':>5}")
    print("-" * 80)

    for m in sorted(mappings, key=lambda x: natural_sort_key(x.target_name, "chr_")):
        print(f"{m.query_name:<15} {m.target_name:<10} {m.coverage:>9.1%} "
              f"{m.plus_strand_length:>12,} {m.minus_strand_length:>12,} "
              f"{'Yes' if m.needs_reverse_complement else 'No':>5}")


def save_mapping_tsv(
    mappings: List[ChromosomeMapping],
    assignments: List[FinalChromosomeAssignment],
    output_path: Path
) -> None:
    """
    Save chromosome mapping summary to TSV file.

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

    print(f"Mapping summary saved to: {output_path}")


def calculate_genome_length(fasta_path: Path) -> int:
    """
    Calculate total genome length from FASTA file.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Total length of all sequences
    """
    sequences = read_fasta(fasta_path)
    return sum(len(seq) for seq in sequences.values())


def write_fasta(
    sequences: Dict[str, str],
    assignments: List[FinalChromosomeAssignment],
    unloc_mappings: List[UnlocMapping],
    output_path: Path,
    output_prefix: str = "SUPER_",
    line_width: int = 60
) -> None:
    """
    Write FASTA file with renamed chromosomes and reverse complement where needed.

    Orientation is taken from a.needs_reverse_complement (chromosomes) and
    unloc.needs_reverse_complement (unlocalized contigs).

    Args:
        sequences: Original sequences dictionary
        assignments: List of FinalChromosomeAssignment (sorted for output)
        unloc_mappings: List of UnlocMapping for unloc contigs
        output_path: Path to output FASTA file
        output_prefix: Prefix for output chromosome names
        line_width: Line width for sequence output (default 60)
    """
    parent_suffix_lookup = {a.original_name: a.new_suffix for a in assignments}
    processed_sequences = set()

    unloc_by_parent: Dict[str, List[UnlocMapping]] = defaultdict(list)
    for unloc in unloc_mappings:
        unloc_by_parent[unloc.parent_chromosome].append(unloc)
    for parent in unloc_by_parent:
        unloc_by_parent[parent].sort(key=lambda x: x.unloc_number)

    with open(output_path, 'w') as f:
        for a in assignments:
            seq = sequences.get(a.original_name, '')
            if not seq:
                print(f"  Warning: No sequence found for {a.original_name}")
                continue

            processed_sequences.add(a.original_name)

            if a.needs_reverse_complement:
                seq = reverse_complement(seq)

            f.write(f">{a.new_name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')

            for unloc in unloc_by_parent.get(a.original_name, []):
                unloc_seq = sequences.get(unloc.contig_name, '')
                if not unloc_seq:
                    print(f"  Warning: No sequence found for {unloc.contig_name}")
                    continue

                processed_sequences.add(unloc.contig_name)

                if unloc.needs_reverse_complement:
                    unloc_seq = reverse_complement(unloc_seq)

                new_parent_suffix = parent_suffix_lookup.get(
                    unloc.parent_chromosome,
                    extract_chromosome_suffix(unloc.parent_chromosome, output_prefix)
                )
                new_name = f"{output_prefix}{new_parent_suffix}_unloc_{unloc.unloc_number}"

                f.write(f">{new_name}\n")
                for i in range(0, len(unloc_seq), line_width):
                    f.write(unloc_seq[i:i+line_width] + '\n')

        for seq_name, seq in sequences.items():
            if seq_name in processed_sequences:
                continue

            f.write(f">{seq_name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')

    print(f"FASTA written to: {output_path}")


def write_chromosome_list(
    assignments: List[FinalChromosomeAssignment],
    unloc_mappings: List[UnlocMapping],
    output_path: Path,
    output_prefix: str = "SUPER_"
) -> None:
    """
    Write chromosome list CSV file.

    Format: name,suffix,yes/no (yes for main chr, no for unloc)

    Args:
        assignments: List of FinalChromosomeAssignment (sorted for output)
        unloc_mappings: List of UnlocMapping for unloc contigs
        output_path: Path to output CSV file
        output_prefix: Prefix for output chromosome names
    """
    parent_suffix_lookup = {a.original_name: a.new_suffix for a in assignments}

    unloc_by_parent_csv: Dict[str, List[UnlocMapping]] = defaultdict(list)
    for unloc in unloc_mappings:
        unloc_by_parent_csv[unloc.parent_chromosome].append(unloc)
    for parent in unloc_by_parent_csv:
        unloc_by_parent_csv[parent].sort(key=lambda x: x.unloc_number)

    with open(output_path, 'w') as f:
        for a in assignments:
            f.write(f"{a.new_name},{a.new_suffix},yes\n")

            for unloc in unloc_by_parent_csv.get(a.original_name, []):
                new_parent_suffix = parent_suffix_lookup.get(
                    unloc.parent_chromosome,
                    extract_chromosome_suffix(unloc.parent_chromosome, output_prefix)
                )
                new_name = f"{output_prefix}{new_parent_suffix}_unloc_{unloc.unloc_number}"
                f.write(f"{new_name},{new_parent_suffix},no\n")

    print(f"Chromosome list written to: {output_path}")
