"""Alignment scatter plots, split into data preparation, rendering and export.

Only :func:`render_alignment_figure` (and the wrapper that uses it) needs
matplotlib, which is imported lazily and never through ``pyplot``, so importing
this module and calling :func:`prepare_plot_data` work without it.
"""
import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from .models import ChromosomeMapping, PAFRecord
from .names import is_unloc_contig

if TYPE_CHECKING:
    from matplotlib.figure import Figure

logger = logging.getLogger(__name__)

PLOT_FORMATS = ("png", "pdf", "svg")


class MissingPlotDependency(ImportError):
    """Raised when matplotlib is required for plotting but is not installed."""


@dataclass(frozen=True)
class PlotData:
    """Plot-ready data for one query->target pair (no matplotlib involved).

    Midpoints are in megabases. ``*_target_mb`` and ``*_query_mb`` are parallel
    tuples, one entry per alignment block.

    Attributes:
        query_name: Query chromosome name.
        target_name: Reference chromosome name.
        target_length: Real length of the target chromosome (bp).
        needs_rc: Orientation decision for the query.
        coverage: Fraction of the query covered by alignments.
        plus_strand_length: Aligned bases on the plus strand.
        minus_strand_length: Aligned bases on the minus strand.
        plus_target_mb: Target midpoints of '+' blocks.
        plus_query_mb: Query midpoints of '+' blocks.
        minus_target_mb: Target midpoints of '-' blocks.
        minus_query_mb: Query midpoints of '-' blocks.
    """

    query_name: str
    target_name: str
    target_length: int
    needs_rc: bool
    coverage: float
    plus_strand_length: int
    minus_strand_length: int
    plus_target_mb: tuple[float, ...]
    plus_query_mb: tuple[float, ...]
    minus_target_mb: tuple[float, ...]
    minus_query_mb: tuple[float, ...]


def prepare_plot_data(
    records: list[PAFRecord],
    mappings: list[ChromosomeMapping],
) -> list[PlotData]:
    """Collect the per-chromosome data needed to draw the alignment plots.

    Unloc queries and mappings without alignment blocks are skipped.

    Args:
        records: Filtered PAF records.
        mappings: Chromosome mappings with orientation decisions.

    Returns:
        One PlotData per plotted mapping, in mapping order.
    """
    blocks_index: dict[str, list[PAFRecord]] = defaultdict(list)
    for r in records:
        blocks_index[f"{r.query_name}->{r.target_name}"].append(r)

    result: list[PlotData] = []
    for mapping in mappings:
        if is_unloc_contig(mapping.query_name):
            continue
        recs = blocks_index.get(f"{mapping.query_name}->{mapping.target_name}", [])
        if not recs:
            continue

        plus_recs = [r for r in recs if r.strand == "+"]
        minus_recs = [r for r in recs if r.strand == "-"]
        result.append(PlotData(
            query_name=mapping.query_name,
            target_name=mapping.target_name,
            # Real target length (not max(target_end)), needed to mirror the right-hand panel.
            target_length=recs[0].target_length,
            needs_rc=mapping.needs_reverse_complement,
            coverage=mapping.coverage,
            plus_strand_length=mapping.plus_strand_length,
            minus_strand_length=mapping.minus_strand_length,
            plus_target_mb=tuple((r.target_start + r.target_end) / 2e6 for r in plus_recs),
            plus_query_mb=tuple((r.query_start + r.query_end) / 2e6 for r in plus_recs),
            minus_target_mb=tuple((r.target_start + r.target_end) / 2e6 for r in minus_recs),
            minus_query_mb=tuple((r.query_start + r.query_end) / 2e6 for r in minus_recs),
        ))
    return result


def _import_figure() -> "type[Figure]":
    """Import matplotlib's Figure lazily, raising MissingPlotDependency if absent."""
    try:
        from matplotlib.figure import Figure
    except ImportError as exc:
        raise MissingPlotDependency(
            "matplotlib not installed -- skipping plots. "
            "Install with: pip install matplotlib"
        ) from exc
    return Figure


def render_alignment_figure(data: PlotData) -> "Figure":
    """Draw the two-panel alignment figure for one chromosome pair.

    Left panel:  raw alignment blocks (target_mid vs query_mid), coloured by strand.
    Right panel: same data after applying the orientation decision (target axis
                 flipped when RC is needed), so a correct decision shows a clean
                 diagonal from bottom-left to top-right.

    Uses ``matplotlib.figure.Figure`` with an Agg canvas (no pyplot, no global state).

    Args:
        data: Prepared plot data.

    Returns:
        The rendered figure.

    Raises:
        MissingPlotDependency: If matplotlib is not installed.
    """
    figure_cls = _import_figure()
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    fig = figure_cls(figsize=(16, 7))
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 2)

    target_len_mb = data.target_length / 1e6

    def scatter(ax, plus_x, minus_x) -> None:  # type: ignore[no-untyped-def]
        if data.plus_query_mb:
            ax.scatter(plus_x, data.plus_query_mb,
                       s=1, alpha=0.4, color="steelblue", label="+ strand")
        if data.minus_query_mb:
            ax.scatter(minus_x, data.minus_query_mb,
                       s=1, alpha=0.4, color="firebrick", label="- strand")
        ax.set_xlabel(f"{data.target_name} position (Mb)")
        ax.set_ylabel(f"{data.query_name} position (Mb)")
        ax.legend(markerscale=8, loc="upper left")

    ax = axes[0]
    scatter(ax, data.plus_target_mb, data.minus_target_mb)
    ax.set_title("All alignment blocks")

    ax2 = axes[1]
    if data.needs_rc:
        scatter(
            ax2,
            [target_len_mb - x for x in data.plus_target_mb],
            [target_len_mb - x for x in data.minus_target_mb],
        )
        ax2.set_title("After RC (target axis flipped)")
    else:
        scatter(ax2, data.plus_target_mb, data.minus_target_mb)
        ax2.set_title("As-is (no RC applied)")

    rc_str = "RC" if data.needs_rc else "no RC"
    fig.suptitle(
        f"{data.query_name}  ->  {data.target_name}   |   "
        f"decision: {rc_str}   coverage: {data.coverage:.1%}   "
        f"+: {data.plus_strand_length:,}  -: {data.minus_strand_length:,}",
        fontsize=11, fontweight="bold",
    )
    fig.tight_layout()
    return fig


def save_figure(fig: "Figure", path: Path, dpi: int = 120) -> None:
    """Save a figure; the file format is taken from the path extension.

    Args:
        fig: Figure to save.
        path: Output file path (.png, .pdf, .svg, ...).
        dpi: Resolution for raster formats.
    """
    fig.savefig(path, dpi=dpi)


def plot_chromosome_alignments(
    records: list[PAFRecord],
    mappings: list[ChromosomeMapping],
    output_dir: Path,
    fmt: str = "png",
) -> list[Path]:
    """Generate scatter plots of PAF alignment blocks for each mapped chromosome.

    Files are written to ``<output_dir>/plots/{query}_vs_{target}.{fmt}``.

    Args:
        records: Filtered PAF records.
        mappings: Chromosome mappings with orientation decisions.
        output_dir: Directory under which the ``plots`` directory is created.
        fmt: Output format: one of png, pdf, svg.

    Returns:
        Paths of the written plot files.

    Raises:
        ValueError: If fmt is not supported.
        MissingPlotDependency: If matplotlib is not installed.
    """
    if fmt not in PLOT_FORMATS:
        raise ValueError(f"Unsupported plot format {fmt!r}; choose from {', '.join(PLOT_FORMATS)}")
    _import_figure()

    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    for data in prepare_plot_data(records, mappings):
        fig = render_alignment_figure(data)
        out_path = plots_dir / f"{data.query_name}_vs_{data.target_name}.{fmt}"
        save_figure(fig, out_path)
        written.append(out_path)

    logger.info("Plots saved to: %s", plots_dir)
    return written
