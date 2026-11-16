"""Tests for report formatting and the logging/stdout split."""
import subprocess
import sys

from rename_and_orient import ChromosomeMapping, FinalChromosomeAssignment
from rename_and_orient.report import format_assignment_table, format_mapping_summary

CLI_CMD = [sys.executable, "-m", "rename_and_orient.cli"]


def _cli(paf_file, fasta_file, out_dir, *extra):
    return subprocess.run(
        CLI_CMD + [
            "--fasta", str(fasta_file), "--paf", str(paf_file),
            "--output-dir", str(out_dir), "--output-prefix", "t", *extra,
        ],
        capture_output=True, text=True,
    )


def test_format_assignment_table():
    a = FinalChromosomeAssignment(
        original_name="SUPER_9", new_name="SUPER_1", new_suffix="1",
        needs_reverse_complement=True, is_sex_chromosome=False,
    )
    text = format_assignment_table([a])
    lines = text.split("\n")
    assert text.startswith("\nFinal chromosome assignments:")
    assert lines[2] == "-" * 70
    assert lines[3] == f"{'Original':<20} {'New Name':<20} {'Suffix':<10} {'RC?':>5} {'Sex?':>5}"
    assert lines[5] == (
        f"{'SUPER_9':<20} {'SUPER_1':<20} {'1':<10} {'Yes':>5} {'No':>5}"
    )


def test_format_mapping_summary():
    m = ChromosomeMapping(
        query_name="SUPER_9", target_name="chr_1", query_length=1000,
        total_alignment_length=900, coverage=0.9,
        plus_strand_length=100, minus_strand_length=800, needs_reverse_complement=True,
    )
    lines = format_mapping_summary([m]).split("\n")
    assert lines[1] == "Chromosome mapping summary:"
    assert lines[2] == "-" * 80
    assert lines[5] == (
        f"{'SUPER_9':<15} {'chr_1':<10} {0.9:>9.1%} {100:>12,} {800:>12,} {'Yes':>5}"
    )


def test_cli_table_on_stdout_progress_on_stderr(paf_file, fasta_file, tmp_path):
    result = _cli(paf_file, fasta_file, tmp_path / "o")
    assert result.returncode == 0, result.stderr
    assert "Final chromosome assignments:" in result.stdout
    assert "Chromosome mapping summary:" in result.stdout
    assert "Reading FASTA file" not in result.stdout
    assert "Reading FASTA file" in result.stderr


def test_cli_quiet_hides_info_keeps_tables(paf_file, fasta_file, tmp_path):
    result = _cli(paf_file, fasta_file, tmp_path / "o", "--quiet")
    assert result.returncode == 0, result.stderr
    assert "Reading FASTA file" not in result.stderr
    assert "Final chromosome assignments:" in result.stdout


def test_cli_log_file_contains_info(paf_file, fasta_file, tmp_path):
    log = tmp_path / "run.log"
    result = _cli(paf_file, fasta_file, tmp_path / "o", "--quiet", "--log-file", str(log))
    assert result.returncode == 0, result.stderr
    text = log.read_text()
    assert " INFO rename_and_orient.cli: Reading FASTA file" in text


def test_library_call_prints_nothing(tmp_path):
    code = (
        "from pathlib import Path\n"
        "import sys\n"
        "from rename_and_orient import read_fasta, write_fasta\n"
        f"p = Path(r'{tmp_path}') / 'a.fa'\n"
        "p.write_text('>SUPER_1\\nACGT\\n')\n"
        "print(len(read_fasta(p)))\n"
        "write_fasta({'SUPER_1': 'ACGT'}, [], [], Path(r'" + str(tmp_path) + "') / 'b.fa')\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "1"
    assert result.stderr == ""
