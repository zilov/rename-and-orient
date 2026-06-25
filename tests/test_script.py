"""Tests for rename-and-orient CLI."""
import subprocess
import sys
import tempfile
from pathlib import Path

HEADER = (
    "query\ttarget\trenamed_to\tquery_length\talignment_length\t"
    "coverage\tplus_strand\tminus_strand\tneeds_reverse_complement\n"
)
TEST_DATA = Path(__file__).parent / "test_data"

# Use 'python -m rename_and_orient.cli' so tests work regardless of PATH.
CLI_CMD = [sys.executable, "-m", "rename_and_orient.cli"]


def _run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True)


def test_script_runs_on_test_data(paf_file, fasta_file):
    """Test that CLI runs successfully on test data."""
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir) / "output"
        output_dir.mkdir()

        result = _run(
            CLI_CMD + [
                "--fasta", str(fasta_file),
                "--paf", str(paf_file),
                "--output-dir", str(output_dir),
                "--output-prefix", "test",
                "--min-coverage", "0.5",
            ]
        )

        assert result.returncode == 0, f"CLI failed with stderr: {result.stderr}"
        assert (output_dir / "test.fa").exists()
        assert (output_dir / "test.chromosome.list.csv").exists()
        assert (output_dir / "test.mapping.tsv").exists()

        content = (output_dir / "test.fa").read_text()
        assert '>' in content
        assert len(content) > 100


def test_mapping_table_mode_produces_output(tmp_path):
    """--mapping-table mode renames FASTA using a pre-built TSV, no PAF needed."""
    fasta = tmp_path / "hap2.fa"
    fasta.write_text(">SUPER_9\nATCGATCGATCG\n>SUPER_3\nGCTAGCTAGCTA\n")

    table = tmp_path / "mapping.tsv"
    table.write_text(
        HEADER
        + "SUPER_9\tchr_1\tSUPER_1\t12\t11\t0.92\t2\t9\tyes\n"
        + "SUPER_3\tchr_2\tSUPER_2\t12\t11\t0.92\t9\t2\tno\n"
    )

    out_dir = tmp_path / "out"
    result = _run(
        CLI_CMD + [
            "--fasta", str(fasta),
            "--mapping-table", str(table),
            "--output-dir", str(out_dir),
            "--output-prefix", "result",
        ]
    )

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert (out_dir / "result.fa").exists()
    assert (out_dir / "result.chromosome.list.csv").exists()

    headers = [
        line.strip()[1:]
        for line in (out_dir / "result.fa").read_text().splitlines()
        if line.startswith(">")
    ]
    assert "SUPER_1" in headers
    assert "SUPER_2" in headers


def test_mapping_table_mode_applies_rc(tmp_path):
    """Chromosome with needs_reverse_complement=yes must be RC'd in output FASTA."""
    seq = "AAAAACCCCC"  # RC = GGGGTTTTT
    fasta = tmp_path / "hap2.fa"
    fasta.write_text(f">SUPER_1\n{seq}\n")

    table = tmp_path / "mapping.tsv"
    table.write_text(HEADER + "SUPER_1\tchr_1\tSUPER_1\t10\t9\t0.9\t1\t8\tyes\n")

    out_dir = tmp_path / "out"
    _run(
        CLI_CMD + [
            "--fasta", str(fasta),
            "--mapping-table", str(table),
            "--output-dir", str(out_dir),
            "--output-prefix", "result",
        ]
    )

    content = (out_dir / "result.fa").read_text()
    assert "GGGGTTTTT" in content, f"Expected RC sequence in output, got:\n{content}"


def test_mapping_table_mode_unloc_follow_parent(tmp_path):
    """Unloc contigs unique to haplotype 2 appear right after their parent in output."""
    fasta = tmp_path / "hap2.fa"
    fasta.write_text(
        ">SUPER_5\nAAAAA\n>SUPER_5_unloc_1\nCCCCC\n>SUPER_3\nGGGGG\n"
    )

    table = tmp_path / "mapping.tsv"
    table.write_text(
        HEADER
        + "SUPER_5\tchr_1\tSUPER_1\t5\t5\t1.0\t5\t0\tno\n"
        + "SUPER_3\tchr_2\tSUPER_2\t5\t5\t1.0\t5\t0\tno\n"
    )

    out_dir = tmp_path / "out"
    _run(
        CLI_CMD + [
            "--fasta", str(fasta),
            "--mapping-table", str(table),
            "--output-dir", str(out_dir),
            "--output-prefix", "result",
        ]
    )

    headers = [
        line.strip()[1:]
        for line in (out_dir / "result.fa").read_text().splitlines()
        if line.startswith(">")
    ]
    assert headers.index("SUPER_1_unloc_1") == headers.index("SUPER_1") + 1, (
        f"Unloc must immediately follow parent. Got order: {headers}"
    )


def test_mapping_table_and_paf_are_mutually_exclusive(tmp_path):
    """Passing both --paf and --mapping-table must produce an error."""
    fasta = tmp_path / "hap.fa"
    fasta.write_text(">SUPER_1\nAAAA\n")
    dummy = tmp_path / "dummy.tsv"
    dummy.write_text("")

    result = _run(
        CLI_CMD + [
            "--fasta", str(fasta),
            "--paf", str(dummy),
            "--mapping-table", str(dummy),
            "--output-dir", str(tmp_path / "out"),
            "--output-prefix", "x",
        ]
    )
    assert result.returncode != 0
