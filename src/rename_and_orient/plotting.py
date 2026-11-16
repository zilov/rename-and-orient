"""Alignment scatter plots, split into data preparation, rendering and export.

Only :func:`render_alignment_figure` (and the wrapper that uses it) needs
matplotlib, which is imported lazily and never through ``pyplot``, so importing
this module and calling :func:`prepare_plot_data` work without it.
"""
import logging
import math
import re
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


OVERVIEW_FILENAME = "alignments_overview.pdf"


def _natural_key(name: str) -> list[object]:
    """Sort key that orders embedded integers numerically (chr_2 before chr_10)."""
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", name)]


def _draw_panel(ax, data: PlotData, after: bool) -> None:  # type: ignore[no-untyped-def]
    """Draw one compact alignment panel (before or after the orientation decision)."""
    target_len_mb = data.target_length / 1e6
    flip = after and data.needs_rc

    def xs(values: tuple[float, ...]) -> list[float]:
        return [target_len_mb - x for x in values] if flip else list(values)

    if data.plus_query_mb:
        ax.scatter(xs(data.plus_target_mb), data.plus_query_mb,
                   s=3.0, alpha=0.5, color="steelblue", linewidths=0, rasterized=True)
    if data.minus_query_mb:
        ax.scatter(xs(data.minus_target_mb), data.minus_query_mb,
                   s=3.0, alpha=0.5, color="firebrick", linewidths=0, rasterized=True)

    if after:
        status = "RC applied" if data.needs_rc else "no RC"
    else:
        status = "needs RC" if data.needs_rc else "ok"
    ax.set_title(f"{data.query_name} -> {data.target_name}  [{status}]",
                 fontsize=7, color="firebrick" if data.needs_rc else "black")
    ax.tick_params(labelsize=6, length=2)
    ax.set_xlim(0, target_len_mb)
    ax.set_ylim(bottom=0)


def render_overview_figure(
    data_list: list[PlotData],
    after: bool,
    ncols: int = 6,
) -> "Figure":
    """Draw one page: a grid with every chromosome pair's alignment blocks.

    Args:
        data_list: Prepared plot data, one entry per panel (drawn in this order).
        after: False for the raw ("before") view; True to apply the orientation
            decision (target axis mirrored where RC is needed).
        ncols: Number of grid columns.

    Returns:
        The rendered figure.

    Raises:
        MissingPlotDependency: If matplotlib is not installed.
    """
    figure_cls = _import_figure()
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.lines import Line2D

    n = max(len(data_list), 1)
    ncols = max(1, min(ncols, n))
    nrows = math.ceil(n / ncols)
    fig = figure_cls(figsize=(4 * ncols, 3.4 * nrows + 0.9))
    FigureCanvasAgg(fig)
    axes = fig.subplots(nrows, ncols, squeeze=False)

    for i, ax in enumerate(axes.flat):
        if i < len(data_list):
            _draw_panel(ax, data_list[i], after)
        else:
            ax.axis("off")

    title = ("After orientation (target axis mirrored where RC is decided)"
             if after else "Before orientation (as aligned)")
    fig.suptitle(title, fontsize=14, fontweight="bold")
    fig.legend(
        handles=[
            Line2D([], [], marker="o", ls="", color="steelblue", markersize=4, label="+ strand"),
            Line2D([], [], marker="o", ls="", color="firebrick", markersize=4, label="- strand"),
        ],
        loc="upper right", fontsize=8,
    )
    fig.supxlabel("target (reference) position, Mb", fontsize=9)
    fig.supylabel("query position, Mb", fontsize=9)
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.97))
    return fig


def save_overview_pdf(data_list: list[PlotData], path: Path, ncols: int = 6) -> None:
    """Write a two-page PDF: page 1 all alignments before, page 2 after orientation.

    Args:
        data_list: Prepared plot data, one entry per panel.
        path: Output PDF path.
        ncols: Number of grid columns.

    Raises:
        MissingPlotDependency: If matplotlib is not installed.
    """
    _import_figure()
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(path) as pdf:
        for after in (False, True):
            fig = render_overview_figure(data_list, after, ncols)
            pdf.savefig(fig)


def plot_alignments_overview(
    records: list[PAFRecord],
    mappings: list[ChromosomeMapping],
    output_dir: Path,
    ncols: int = 6,
) -> Path | None:
    """Write ``<output_dir>/plots/alignments_overview.pdf`` (before / after grids).

    Panels are ordered by target (reference) chromosome name, naturally sorted.

    Args:
        records: Filtered PAF records.
        mappings: Chromosome mappings with orientation decisions.
        output_dir: Directory under which the ``plots`` directory is created.
        ncols: Number of grid columns.

    Returns:
        Path of the written PDF, or None if there is nothing to plot.

    Raises:
        MissingPlotDependency: If matplotlib is not installed.
    """
    _import_figure()
    data_list = sorted(
        prepare_plot_data(records, mappings),
        key=lambda d: (_natural_key(d.target_name), _natural_key(d.query_name)),
    )
    if not data_list:
        logger.warning("No alignments to plot; overview PDF not written")
        return None

    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    out_path = plots_dir / OVERVIEW_FILENAME
    save_overview_pdf(data_list, out_path, ncols)
    logger.info("Alignment overview saved to: %s", out_path)
    return out_path
