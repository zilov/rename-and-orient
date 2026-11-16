"""Tests for the indexed, streaming FASTA reader/writer."""
import gzip
import logging
import tempfile
import tracemalloc
from pathlib import Path

import pytest

import rename_and_orient.fasta as fasta_mod
from rename_and_orient.fasta import (
    FastaIndex,
    calculate_genome_length,
    index_fasta,
    read_fasta,
    write_fasta,
)
from rename_and_orient.models import FinalChromosomeAssignment, UnlocMapping


def _old_read_fasta(path: Path) -> dict[str, str]:
    """Reference implementation: the pre-index reader (text mode, strip every line)."""
    seqs: dict[str, str] = {}
    name = None
    parts: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def _make(tmp_path: Path, text: str, name: str = "x.fa") -> Path:
    path = tmp_path / name
    path.write_bytes(text.encode())
    return path


class TestIndexFasta:
    """index_fasta / FastaIndex behaviour."""

    def test_basic(self, tmp_path):
        path = _make(tmp_path, ">a\nACGT\nAC\n>b\nGGG\n")
        with index_fasta(path) as idx:
            assert idx.names == ["a", "b"]
            assert len(idx) == 2
            assert "a" in idx and "c" not in idx
            assert idx.length("a") == 6 and idx.length("b") == 3
            assert idx.total_length == 9
            assert idx.read_sequence("a") == "ACGTAC"
            assert idx.read_sequence("b") == "GGG"
            assert idx.get("zzz", "") == ""
            assert idx["b"] == "GGG"
            with pytest.raises(KeyError):
                idx.read_sequence("zzz")

    def test_multi_width_and_last_short_line(self, tmp_path):
        text = ">a\n" + "A" * 70 + "\n" + "C" * 70 + "\nGG\n>b\n" + "T" * 60 + "\n" + "T" * 5 + "\n"
        path = _make(tmp_path, text)
        with index_fasta(path) as idx:
            assert idx.read_sequence("a") == "A" * 70 + "C" * 70 + "GG"
            assert idx.read_sequence("b") == "T" * 65

    def test_crlf(self, tmp_path):
        path = _make(tmp_path, ">a desc\r\nACGT\r\nAC\r\n\r\n>b\r\nGG\r\n")
        with index_fasta(path) as idx:
            assert idx.names == ["a", "b"]
            assert idx.read_sequence("a") == "ACGTAC"
            assert idx.length("a") == 6
            assert idx.read_sequence("b") == "GG"

    def test_header_description(self, tmp_path):
        path = _make(tmp_path, ">chr1 length=5 species=foo\nACGTA\n>chr2\tTab desc\nCC\n")
        with index_fasta(path) as idx:
            assert idx.names == ["chr1", "chr2"]

    def test_blank_lines(self, tmp_path):
        path = _make(tmp_path, "\n>a\n\nAC\n\nGT\n\n\n>b\n\nT\n")
        with index_fasta(path) as idx:
            assert idx.read_sequence("a") == "ACGT"
            assert idx.length("a") == 4
            assert idx.read_sequence("b") == "T"

    def test_empty_record(self, tmp_path):
        path = _make(tmp_path, ">a\n>b\nAC\n>c\n")
        with index_fasta(path) as idx:
            assert idx.names == ["a", "b", "c"]
            assert idx.read_sequence("a") == "" and idx.length("a") == 0
            assert idx.read_sequence("c") == ""
            assert idx.read_sequence("b") == "AC"

    def test_whitespace_edges_match_old_reader(self, tmp_path):
        text = ">a\n  ACGT \t\n AC\nG T\n   \n>b\n\tTT\n>c\nGG\n"
        path = _make(tmp_path, text)
        expected = _old_read_fasta(path)
        with index_fasta(path) as idx:
            assert {n: idx.read_sequence(n) for n in idx.names} == expected
            assert idx.read_sequence("a") == "ACGTACG T"  # interior space is kept, as before
            assert idx.length("a") == len(expected["a"])
            assert idx.total_length == sum(len(s) for s in expected.values())

    def test_matches_old_reader_on_mixed_file(self, tmp_path):
        text = (
            ">a desc\r\nacgtNN\r\n\r\nAC\r\n>e\r\n>b\n GG \nT\n\n>c\n"
            + "A" * 61 + "\n" + "C" * 3 + "\n"
        )
        path = _make(tmp_path, text)
        with index_fasta(path) as idx:
            assert {n: idx.read_sequence(n) for n in idx.names} == _old_read_fasta(path)
            assert read_fasta(path) == _old_read_fasta(path)

    def test_duplicate_names_raise(self, tmp_path):
        path = _make(tmp_path, ">a\nAC\n>b\nGG\n>a extra\nTT\n")
        with pytest.raises(ValueError, match="duplicate sequence name 'a'"):
            index_fasta(path)
        with pytest.raises(ValueError, match="duplicate"):
            read_fasta(path)

    def test_gzip_input_uses_temp_file_and_cleans_up(self, tmp_path, monkeypatch):
        tmpdir = tmp_path / "tmp"
        tmpdir.mkdir()
        monkeypatch.setattr(tempfile, "tempdir", str(tmpdir))
        gz = tmp_path / "x.fa.gz"
        with gzip.open(gz, "wt") as f:
            f.write(">a\nACGT\nAC\n>b\nGG\n")
        with index_fasta(gz) as idx:
            assert len(list(tmpdir.iterdir())) == 1  # decompressed copy exists
            assert idx.read_sequence("a") == "ACGTAC"
            assert idx.total_length == 8
        assert list(tmpdir.iterdir()) == []

    def test_gzip_temp_file_removed_on_error(self, tmp_path, monkeypatch):
        tmpdir = tmp_path / "tmp"
        tmpdir.mkdir()
        monkeypatch.setattr(tempfile, "tempdir", str(tmpdir))
        gz = tmp_path / "dup.fa.gz"
        with gzip.open(gz, "wt") as f:
            f.write(">a\nAC\n>a\nGG\n")
        with pytest.raises(ValueError):
            index_fasta(gz)
        assert list(tmpdir.iterdir()) == []

    def test_plain_file_has_no_temp_copy(self, tmp_path, monkeypatch):
        tmpdir = tmp_path / "tmp"
        tmpdir.mkdir()
        monkeypatch.setattr(tempfile, "tempdir", str(tmpdir))
        path = _make(tmp_path, ">a\nAC\n")
        with index_fasta(path) as idx:
            assert idx.read_sequence("a") == "AC"
            assert list(tmpdir.iterdir()) == []


def _synthetic_fasta(tmp_path: Path) -> Path:
    """FASTA with chromosomes, unlocs, shrapnel, an empty record, mixed widths and CRLF."""
    recs = [
        ("SUPER_1", 701, 70), ("SUPER_2", 500, 60), ("SUPER_2_unloc_1", 130, 50),
        ("SUPER_2_unloc_2", 77, 60), ("scaffold_5", 90, 60), ("empty", 0, 60),
        ("SUPER_3", 333, 80),
    ]
    out = []
    for i, (name, n, w) in enumerate(recs):
        seq = ("ACGTNacgtnRY"[(i % 3):] * (n // 8 + 2))[:n]
        out.append(f">{name} d\r\n" if i == 4 else f">{name}\n")
        out.extend(seq[j:j + w] + "\n" for j in range(0, n, w))
    return _make(tmp_path, "".join(out), "syn.fa")


class TestWriteFastaSources:
    """write_fasta output must be identical for dict and FastaIndex sources."""

    @pytest.mark.parametrize("width", [60, 70, 1000])
    def test_dict_and_index_output_identical(self, tmp_path, width, caplog):
        path = _synthetic_fasta(tmp_path)
        asg = [
            FinalChromosomeAssignment("SUPER_2", "SUPER_1", "1", True, False),
            FinalChromosomeAssignment("SUPER_1", "SUPER_2", "2", False, False),
            FinalChromosomeAssignment("missing", "SUPER_9", "9", False, False),
            FinalChromosomeAssignment("empty", "SUPER_8", "8", True, False),
        ]
        unl = [
            UnlocMapping("SUPER_2_unloc_2", "SUPER_2", 2, True),
            UnlocMapping("SUPER_2_unloc_1", "SUPER_2", 1, False),
        ]
        out_dict, out_idx = tmp_path / "d.fa", tmp_path / "i.fa"
        with caplog.at_level(logging.WARNING):
            write_fasta(_old_read_fasta(path), asg, unl, out_dict, "SUPER_", width)
            n_dict = len(caplog.records)
            with index_fasta(path) as idx:
                write_fasta(idx, asg, unl, out_idx, "SUPER_", width)
        assert out_dict.read_bytes() == out_idx.read_bytes()
        assert len(caplog.records) == 2 * n_dict > 0  # same warnings for both sources
        text = out_idx.read_text()
        assert ">SUPER_1\n" in text and "scaffold_5" in text and ">empty\n" in text


class TestStreamingLength:
    """calculate_genome_length streams and agrees with the index."""

    def test_plain_and_gzip(self, tmp_path, monkeypatch):
        path = _synthetic_fasta(tmp_path)
        gz = tmp_path / "syn.fa.gz"
        gz.write_bytes(gzip.compress(path.read_bytes()))
        monkeypatch.setattr(fasta_mod, "read_fasta", None)  # must not be used
        with index_fasta(path) as idx:
            expected = idx.total_length
        assert expected == 701 + 500 + 130 + 77 + 90 + 0 + 333
        assert calculate_genome_length(path) == expected
        assert calculate_genome_length(gz) == expected

    def test_ignores_text_before_first_header(self, tmp_path):
        path = _make(tmp_path, "ACGT\n>a\nAC\n")
        assert calculate_genome_length(path) == 2


class TestIndexIsLazy:
    """Reading through the index must not materialise all records."""

    def test_write_reads_each_record_once(self, tmp_path, monkeypatch):
        path = _synthetic_fasta(tmp_path)
        calls: list[str] = []
        orig = FastaIndex.read_sequence

        def spy(self, name):
            calls.append(name)
            return orig(self, name)

        monkeypatch.setattr(FastaIndex, "read_sequence", spy)
        with index_fasta(path) as idx:
            assert calls == []  # indexing reads no sequences
            write_fasta(idx, [], [], tmp_path / "o.fa")
        assert sorted(calls) == sorted(["SUPER_1", "SUPER_2", "SUPER_2_unloc_1", "SUPER_2_unloc_2",
                                        "scaffold_5", "empty", "SUPER_3"])

    def test_peak_memory_well_below_file_size(self, tmp_path):
        n_rec, rec_len = 40, 250_000
        with open(tmp_path / "big.fa", "w") as f:
            line = "ACGT" * 15 + "\n"
            for i in range(n_rec):
                f.write(f">c{i}\n")
                f.write(line * (rec_len // 60))
        total = n_rec * (rec_len // 60) * 60
        out = tmp_path / "out.fa"
        asg = [FinalChromosomeAssignment("c0", "SUPER_1", "1", True, False)]
        tracemalloc.start()
        try:
            with index_fasta(tmp_path / "big.fa") as idx:
                write_fasta(idx, asg, [], out)
            peak = tracemalloc.get_traced_memory()[1]
        finally:
            tracemalloc.stop()
        assert out.stat().st_size > total
        assert peak < total / 4, f"peak {peak} vs total {total}"
