"""Orientation decision (reverse-complement or not) from alignment blocks."""
from collections.abc import Sequence

from .models import PAFRecord


def _pearson_r(
    xs: Sequence[float], ys: Sequence[float], ws: Sequence[float] | None = None
) -> float:
    """Compute the weighted Pearson correlation coefficient.

    Weights default to 1 (plain Pearson). Returns 0.0 if fewer than 2 points
    or zero variance in either variable.

    Args:
        xs: First variable.
        ys: Second variable (same length as xs).
        ws: Optional per-point weights (same length as xs).

    Returns:
        Correlation coefficient in [-1, 1], or 0.0 if undefined.
    """
    n = len(xs)
    if n < 2:
        return 0.0
    if ws is None:
        ws = [1.0] * n
    total_w = sum(ws)
    if total_w == 0:
        return 0.0
    mx = sum(w * x for w, x in zip(ws, xs)) / total_w
    my = sum(w * y for w, y in zip(ws, ys)) / total_w
    num = sum(w * (x - mx) * (y - my) for w, x, y in zip(ws, xs, ys))
    sx = sum(w * (x - mx) ** 2 for w, x in zip(ws, xs)) ** 0.5
    sy = sum(w * (y - my) ** 2 for w, y in zip(ws, ys)) ** 0.5
    if sx == 0 or sy == 0:
        return 0.0
    return num / (sx * sy)


def needs_reverse_complement_by_correlation(
    records: list[PAFRecord],
    target_length: int,
) -> bool:
    """Decide whether a query chromosome needs reverse-complement.

    Uses Pearson correlation between query midpoints and target midpoints.

    For '+' blocks: r = corr(query_mid, target_mid).
    For '-' blocks: r = corr(query_mid, target_length - target_mid).

    Each block is weighted by its query length, so that many short scattered
    repeat hits cannot outvote a few long syntenic blocks.

    The strand whose r is higher wins. This correctly ignores large structural
    inversions (their blocks are scattered -> low r) while syntenic blocks
    form a clean diagonal -> high r.

    Args:
        records:       PAFRecord objects for one query->best_target pair.
        target_length: Length of the target chromosome.

    Returns:
        True  if reverse complement is needed,
        False if the chromosome is already correctly oriented.
    """
    plus_recs = [r for r in records if r.strand == '+']
    minus_recs = [r for r in records if r.strand == '-']

    r_plus = 0.0
    if len(plus_recs) >= 2:
        xs = [(r.query_start + r.query_end) / 2 for r in plus_recs]
        ys = [(r.target_start + r.target_end) / 2 for r in plus_recs]
        ws = [r.query_end - r.query_start for r in plus_recs]
        r_plus = _pearson_r(xs, ys, ws)

    r_minus = 0.0
    if len(minus_recs) >= 2:
        xs = [(r.query_start + r.query_end) / 2 for r in minus_recs]
        ys = [target_length - (r.target_start + r.target_end) / 2 for r in minus_recs]
        ws = [r.query_end - r.query_start for r in minus_recs]
        r_minus = _pearson_r(xs, ys, ws)

    # If correlation is uninformative for both strands (< 2 blocks each),
    # fall back to naive length comparison.
    if r_plus == 0.0 and r_minus == 0.0:
        plus_len = sum(r.query_end - r.query_start for r in plus_recs)
        minus_len = sum(r.query_end - r.query_start for r in minus_recs)
        return minus_len > plus_len

    return r_minus > r_plus
