"""Tests for the module split, the deprecated io shim and the plotting refactor."""
import importlib
import subprocess
import sys

import pytest

from rename_and_orient import ChromosomeMapping, PAFRecord
from rename_and_orient.plotting import (
    MissingPlotDependency,
    PlotData,
    plot_chromosome_alignments,
    prepare_plot_data,
    render_alignment_figure,
    save_figure,
)

TARGET_LEN = 50_000_000


def _rec(query, strand, qs, qe, ts, te, target="chr_1"):
    return PAFRecord(query, 10_000_000, qs, qe, strand, target, TARGET_LEN, ts, te,
                     qe - qs, qe - qs, 60)


def _mapping(query, target="chr_1", rc=False):
    return ChromosomeMapping(query, 10_000_000, target, 4_000_000, 0.4, 2_000_000, 2_000_000, rc)


@pytest.fixture
def synthetic():
    records = [
        _rec("SUPER_1", "+", 0, 2_000_000, 1_000_000, 3_000_000),
        _rec("SUPER_1", "-", 4_000_000, 6_000_000, 10_000_000, 12_000_000),
        _rec("SUPER_1_unloc_1", "+", 0, 100_000, 0, 100_000),
        _rec("SUPER_2", "+", 0, 1_000_000, 0, 1_000_000, target="chr_2"),
    ]
    mappings = [
        _mapping("SUPER_1", rc=True),
        _mapping("SUPER_1_unloc_1"),
        _mapping("SUPER_3", target="chr_3"),  # no blocks -> skipped
        _mapping("SUPER_2", target="chr_2"),
    ]
    return records, mappings


def test_io_shim_warns_and_exposes_old_names():
    sys.modules.pop("rename_and_orient.io", None)
    with pytest.warns(DeprecationWarning, match="rename_and_orient.io is deprecated"):
        io = importlib.import_module("rename_and_orient.io")
    for name in (
        "read_fasta", "write_fasta", "calculate_genome_length", "reverse_complement",
        "parse_paf", "natural_sort_key", "print_mapping_summary", "save_mapping_tsv",
        "write_chromosome_list", "_unloc_output_name", "_RC_TABLE",
    ):
        assert hasattr(io, name)


def test_public_api_and_cli_emit_no_deprecation_warning(paf_file, fasta_file, tmp_path):
    code = (
        "import warnings; warnings.simplefilter('error', DeprecationWarning);"
        "import rename_and_orient, rename_and_orient.cli, rename_and_orient.alignment"
    )
    assert subprocess.run([sys.executable, "-c", code]).returncode == 0
    result = subprocess.run(
        [sys.executable, "-W", "error::DeprecationWarning", "-m", "rename_and_orient.cli",
         "--fasta", str(fasta_file), "--paf", str(paf_file), "--output-dir", str(tmp_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr


def test_alignment_reexports_moved_names():
    from rename_and_orient import alignment
    for name in (
        "needs_reverse_complement_by_correlation", "filter_paf_records",
        "detect_reference_prefix", "validate_paf_fasta_consistency",
        "plot_chromosome_alignments",
    ):
        assert callable(getattr(alignment, name))


class TestPreparePlotData:
    def test_skips_unloc_and_mappings_without_blocks(self, synthetic):
        data = prepare_plot_data(*synthetic)
        assert [d.query_name for d in data] == ["SUPER_1", "SUPER_2"]

    def test_midpoints_in_mb_and_strand_split(self, synthetic):
        d = prepare_plot_data(*synthetic)[0]
        assert isinstance(d, PlotData)
        assert d.plus_target_mb == (2.0,)
        assert d.plus_query_mb == (1.0,)
        assert d.minus_target_mb == (11.0,)
        assert d.minus_query_mb == (5.0,)

    def test_target_len_comes_from_records(self, synthetic):
        d = prepare_plot_data(*synthetic)[0]
        assert d.target_length == TARGET_LEN
        assert d.needs_rc is True
        assert d.target_name == "chr_1"

    def test_is_frozen(self, synthetic):
        d = prepare_plot_data(*synthetic)[0]
        with pytest.raises(AttributeError):
            d.query_name = "x"  # type: ignore[misc]


class TestRendering:
    @pytest.mark.parametrize("ext", ["png", "svg", "pdf"])
    def test_render_and_save(self, synthetic, tmp_path, ext):
        pytest.importorskip("matplotlib")
        fig = render_alignment_figure(prepare_plot_data(*synthetic)[0])
        out = tmp_path / f"fig.{ext}"
        save_figure(fig, out)
        assert out.stat().st_size > 0

    @pytest.mark.parametrize("fmt", ["png", "svg", "pdf"])
    def test_wrapper_file_names(self, synthetic, tmp_path, fmt):
        pytest.importorskip("matplotlib")
        paths = plot_chromosome_alignments(*synthetic, tmp_path, fmt)
        assert [p.name for p in paths] == [
            f"SUPER_1_vs_chr_1.{fmt}", f"SUPER_2_vs_chr_2.{fmt}",
        ]
        assert all(p.parent == tmp_path / "plots" and p.exists() for p in paths)

    def test_wrapper_rejects_unknown_format(self, synthetic, tmp_path):
        with pytest.raises(ValueError, match="Unsupported plot format"):
            plot_chromosome_alignments(*synthetic, tmp_path, "gif")


def test_missing_matplotlib_raises_dedicated_error(synthetic, tmp_path, monkeypatch):
    for name in [m for m in sys.modules if m == "matplotlib" or m.startswith("matplotlib.")]:
        monkeypatch.setitem(sys.modules, name, None)
    monkeypatch.setitem(sys.modules, "matplotlib", None)
    assert issubclass(MissingPlotDependency, ImportError)
    with pytest.raises(MissingPlotDependency):
        render_alignment_figure(prepare_plot_data(*synthetic)[0])
    with pytest.raises(MissingPlotDependency):
        plot_chromosome_alignments(*synthetic, tmp_path)
    assert not (tmp_path / "plots").exists()


def test_cli_plot_alignments(paf_file, fasta_file, tmp_path):
    pytest.importorskip("matplotlib")
    result = subprocess.run(
        [sys.executable, "-m", "rename_and_orient.cli", "--fasta", str(fasta_file),
         "--paf", str(paf_file), "--output-dir", str(tmp_path), "-P", "--plot-format", "svg"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert list((tmp_path / "plots").glob("*.svg"))


def test_cli_warns_without_matplotlib(paf_file, fasta_file, tmp_path):
    """Simulate a missing matplotlib in a subprocess; the run continues with a warning."""
    code = (
        "import sys; sys.modules['matplotlib'] = None;"
        "f, p, d = sys.argv[1:4];"
        "sys.argv = ['x', '--fasta', f, '--paf', p, '--output-dir', d, '-P'];"
        "from rename_and_orient.cli import main; main()"
    )
    result = subprocess.run(
        [sys.executable, "-c", code, str(fasta_file), str(paf_file), str(tmp_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "matplotlib not installed" in result.stderr
