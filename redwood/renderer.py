"""Modern redwood circular plotting backend."""

from __future__ import annotations

import json
import re
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/redwood-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.colors import to_rgba
from matplotlib.patches import Polygon, Wedge
import numpy as np
import pysam

from .functions import print_images
from .multipass import (
    MULTIPASS_COLORS,
    classify_circular_reads,
    pack_circular_reads,
)


PLOT_LIMIT = 1.30

FEATURE_COLORS = {
    "gene": "#2f8f46",
    "CDS": "#2f8f46",
    "rRNA": "#c23b3b",
    "tRNA": "#d870a2",
}

from .numts import CLASS_COLORS  # read-class colors shared with the NUMT module

# IGV-style read marks: mismatched read base, insertion, deletion.
MARK_COLORS = {"A": "#009900", "C": "#0000ff", "G": "#d17105", "T": "#ff0000", "I": "#800080", "D": "#000000"}
VARIANT_RING_COLORS = {"mismatch": "#d62728", "deletion": "#000000", "insertion": "#800080", "minor": "#f28e2b"}

# Regular read arcs: nominal radial thickness, plus per-CIGAR-op widths so a
# read visibly shows insertions (fat bulge) and deletions (thin near-gap).
READ_ARC_WIDTH = 0.0058
CIGAR_OP_WIDTH = {
    1: 0.0082,   # I  insertion — fat bulge (kept under the 0.0087 rung spacing)
    2: 0.0020,   # D  deletion — thin near-gap
    3: 0.0020,   # N  reference skip — thin near-gap
}

# Multi-pass spirals are capped at this fraction of the read band's radius;
# the remainder of the band goes to regular single-pass reads.
MULTIPASS_RADIUS_FRACTION = 0.40

# Redwood-wood gradient for regular read arcs — dark heartwood at the read's
# genomic start, through cinnamon, to pale sapwood at its end. The color
# sweep makes read direction visible.
REDWOOD_GRADIENT = ["#3a1d10", "#9c5a2c", "#d3a878"]

# Tree-anatomy theming: the RNA-seq depth ring is recolored as bark, and the
# AT-content ring uses a warm YlOrBr colormap — dynamic-range normalized in
# add_at_track so subtle AT% differences stay readable.
BARK_COLOR = "#6e4a33"
BARK_COLOR_ALT = "#9a6a44"
AT_COLORMAP = plt.cm.YlOrBr
# The AT ring's dynamic range maps into this sub-range of the colormap: the
# floor keeps the lightest (low-AT) values from washing out toward white, and
# using less than the full palette keeps the contrast from being harsh.
AT_RANGE = (0.38, 0.92)


def theta(pos: int, length: int) -> float:
    return 90 - ((pos % length) / length * 360)


def polar_xy(radius: float, angle: float) -> tuple[float, float]:
    radians = np.deg2rad(angle)
    return radius * np.cos(radians), radius * np.sin(radians)


def _add_patch_fast(ax, patch, sink: list | None = None) -> None:
    """Add a patch without the per-patch data-limit update of ``Axes.add_patch``.

    Redwood fixes the axis limits itself, and ``add_patch`` recomputes the Bezier extrema of every wedge to update the data
    limits, which dominated render time (a map has ~10^4 read wedges). With ``sink`` the patch is collected for a
    ``PatchCollection`` instead."""
    if sink is not None:
        sink.append(patch)
    else:
        ax.add_artist(patch)


def add_arc(ax, start: int, stop: int, length: int, radius: float, width: float, sink: list | None = None, **kwargs) -> None:
    if stop <= start:
        return
    if stop - start >= length:
        _add_patch_fast(ax, Wedge((0, 0), radius, 0, 360, width=width, **kwargs), sink)
        return
    if stop > length:
        add_arc(ax, start, length, length, radius, width, sink=sink, **kwargs)
        add_arc(ax, 0, stop - length, length, radius, width, sink=sink, **kwargs)
        return
    _add_patch_fast(ax, Wedge((0, 0), radius, theta(stop, length), theta(start, length), width=width, **kwargs), sink)


def add_directional_feature(
    ax,
    start: int,
    stop: int,
    length: int,
    radius: float,
    width: float,
    strand: str,
    color: str,
    alpha: float = 0.95,
) -> None:
    span = stop - start
    if span <= 0:
        return
    if span < 220:
        add_arc(ax, start, stop, length, radius, width, color=color, alpha=alpha, linewidth=0)
        return
    arrow_bp = int(min(max(length * 0.0045, 70), span * 0.28, 240))
    inner = radius - width
    middle = radius - (width / 2)
    if strand == "-":
        add_arc(ax, start + arrow_bp, stop, length, radius, width, color=color, alpha=alpha, linewidth=0)
        base = theta(start + arrow_bp, length)
        tip = theta(start, length)
    else:
        add_arc(ax, start, stop - arrow_bp, length, radius, width, color=color, alpha=alpha, linewidth=0)
        base = theta(stop - arrow_bp, length)
        tip = theta(stop, length)
    ax.add_patch(
        Polygon(
            [polar_xy(inner, base), polar_xy(radius, base), polar_xy(middle, tip)],
            closed=True,
            facecolor=color,
            edgecolor="none",
            alpha=alpha,
        )
    )


def _cigar_width_profile(
    cigar: list[tuple[int, int]], min_indel: int, base_lw: float
) -> np.ndarray:
    """Per-reference-base line widths for a multi-pass spiral: ``base_lw`` for
    matches, thin for deletions, fat for insertions (centered, back-offset by
    half their length). Indels shorter than ``min_indel`` bp keep the match
    width. Insertions are overlaid last so a following match cannot clobber
    them."""
    thin = base_lw * 0.34
    fat = base_lw * 1.5
    span = sum(ln for op, ln in cigar if op in (0, 2, 3, 7, 8)) or 1
    widths = np.full(span, base_lw, dtype=float)
    insertions = []
    ref = 0
    for op, oplen in cigar:
        if op == 1:                          # insertion — overlaid in pass 2
            insertions.append((ref, oplen))
            continue
        if op in (4, 5, 6):                  # S / H / P — no reference span
            continue
        if op in (2, 3) and oplen >= min_indel:
            widths[ref:ref + oplen] = thin
        ref += oplen
    for ref_off, oplen in insertions:
        if oplen >= min_indel:
            lo = max(0, ref_off - oplen // 2)
            hi = min(span, ref_off + oplen - oplen // 2)
            widths[lo:hi] = fat
    return widths


def add_spiral_read(
    ax,
    start_pos: int,
    passes: float,
    length: int,
    r_outer: float,
    rung_width: float,
    wrap_ramp: float,
    color: str,
    cigar: list[tuple[int, int]] | None = None,
    min_indel: int = 10,
    base_linewidth: float | None = None,
    points_per_turn: int = 540,
    marks: list[tuple[int, str]] | None = None,
    mark_filter: set[int] | None = None,
    mark_indels: bool = True,
) -> None:
    """Draw a multi-pass (rolling-circle) read as an inward spiral.

    The read traverses the circle ``passes`` times. Each completed turn drops
    one rung; the radius is flat for the first ``1 - wrap_ramp`` of a turn and
    then ramps diagonally down over the final ``wrap_ramp`` fraction, so the
    step into the next rung is visible. When a CIGAR is supplied the spiral's
    line width varies along its length — fat at insertions, thin at deletions —
    so multi-pass reads carry the same indel detail as the regular reads. ``base_linewidth`` defaults to 60 % of the
    rung spacing in points (capped at 0.9 pt), so spirals stay distinct lines at any figure size.
    """
    if base_linewidth is None:
        # never wider than ~60 % of the rung spacing, so neighboring turns and reads stay separable at any figure size
        base_linewidth = max(0.25, min(0.9, 0.6 * rung_width / units_per_point(ax)))
    n_steps = max(2, int(points_per_turn * passes) + 1)
    a0 = theta(start_pos, length)
    xs, ys = [], []
    for i in range(n_steps):
        t = passes * i / (n_steps - 1)          # turns elapsed, 0..passes
        k = int(max(0.0, t - 1e-9))             # current turn index
        f = t - k                               # fraction within the turn
        if f <= 1.0 - wrap_ramp:
            r = r_outer - rung_width * k
        else:
            ramp = (f - (1.0 - wrap_ramp)) / wrap_ramp
            r = r_outer - rung_width * (k + ramp)
        x, y = polar_xy(r, a0 - t * 360.0)
        xs.append(x)
        ys.append(y)
    if cigar:
        widths = _cigar_width_profile(cigar, min_indel, base_linewidth)
        points = np.column_stack([xs, ys])
        segments = np.stack([points[:-1], points[1:]], axis=1)
        seg_lw = []
        for i in range(len(segments)):
            t_mid = passes * (i + 0.5) / (n_steps - 1)
            offset = min(int(t_mid * length), len(widths) - 1)
            seg_lw.append(float(widths[offset]))
        ax.add_collection(
            LineCollection(
                segments,
                colors=color,
                linewidths=seg_lw,
                capstyle="round",
                joinstyle="round",
                zorder=3,
            )
        )
    else:
        ax.plot(
            xs,
            ys,
            color=color,
            linewidth=base_linewidth,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=3,
        )
    if marks:
        span_bp = passes * length
        for offset, kind in marks:
            if offset < 0 or offset >= span_bp:
                continue
            if kind in ("I", "D"):
                if not mark_indels:
                    continue
            elif mark_filter is not None and ((start_pos + offset) % length) not in mark_filter:
                continue
            t = offset / length
            k = int(max(0.0, t - 1e-9))
            f = t - k
            r = r_outer - rung_width * (k if f <= 1.0 - wrap_ramp else k + (f - (1.0 - wrap_ramp)) / wrap_ramp)
            x, y = polar_xy(r, a0 - t * 360.0)
            ax.plot([x], [y], marker="o", ms=base_linewidth * 1.4, mew=0, color=MARK_COLORS[kind], zorder=4)


def add_read_marks(
    ax,
    start: int,
    marks: list[tuple[int, str]],
    length: int,
    radius: float,
    width: float = READ_ARC_WIDTH,
    mark_filter: set[int] | None = None,
    mark_indels: bool = True,
) -> int:
    """Overlay IGV-style marks on a regular read arc drawn at ``radius`` (outer edge). ``marks`` are
    ``(offset, kind)`` from :func:`redwood.multipass.read_marks`; mismatch marks are restricted to
    ``mark_filter`` positions when given (0-based, folded). Returns the number of marks drawn."""
    half = max(2.0, length / 1400.0)      # keep a mark visible at any genome length
    centerline = radius - width / 2
    drawn = 0
    for offset, kind in marks:
        if offset < 0 or offset >= length:
            continue
        pos = (start + offset) % length
        if kind in ("I", "D"):
            if not mark_indels:
                continue
        elif mark_filter is not None and pos not in mark_filter:
            continue
        add_arc(ax, int(pos - half), int(pos + half) + 1, length, centerline + width / 2, width,
                color=MARK_COLORS[kind], alpha=1.0, linewidth=0, zorder=4)
        drawn += 1
    return drawn


def add_numt_ring(ax, loci: list[dict], length: int, r_in: float, r_out: float, lanes: int = 3) -> int:
    """Mitogenome intervals that exist as NUMTs in the nuclear genome, one arc per locus interval, colored by identity
    (viridis, 75-100 %), stacked into ``lanes`` lanes (largest loci outermost)."""
    from .numts import _identity_color
    lane_h = (r_out - r_in) / lanes; drawn = 0
    ends = [-1] * lanes
    for r in sorted(loci, key=lambda x: -x["mito_bp"]):
        for s, e in r["intervals"]:
            lane = next((i for i in range(lanes) if ends[i] <= s), None)
            if lane is None: lane = min(range(lanes), key=lambda i: ends[i])
            ends[lane] = e
            add_arc(ax, s, e, length, r_out - lane * lane_h, lane_h * 0.85, color=_identity_color(r["identity"]), alpha=0.95, linewidth=0, zorder=2)
            drawn += 1
    return drawn


def add_variant_ring(ax, rows: list[dict[str, object]], length: int, r_in: float = 0.958, r_out: float = 0.998) -> int:
    """Per-column disagreement ring: one bar per flagged column, height = minor-allele fraction
    (full height for columns whose majority differs from the reference), colored by event."""
    half = max(2.0, length / 1400.0)
    drawn = 0
    ax.add_patch(plt.Circle((0, 0), r_in, fill=False, lw=0.3, color="#c8ced8", zorder=1))
    for row in rows:
        events = str(row.get("events", "")).split(",")
        if "mismatch" in events or "deletion" in events:
            kind, frac = ("deletion" if "deletion" in events else "mismatch"), 1.0
        elif "insertion" in events:
            kind, frac = "insertion", max(float(row.get("ins_frac", 0.0)), float(row.get("minor_frac", 0.0)))
        elif "minor" in events:
            kind, frac = "minor", float(row.get("minor_frac", 0.0))
        else:
            continue
        pos = int(row["pos"]) - 1
        h = (r_out - r_in) * max(0.15, min(1.0, frac))
        add_arc(ax, int(pos - half), int(pos + half) + 1, length, r_in + h, h,
                color=VARIANT_RING_COLORS[kind], alpha=0.95, linewidth=0, zorder=2)
        drawn += 1
    return drawn


def load_variant_table(path: Path) -> list[dict[str, object]]:
    rows = []
    lines = Path(path).read_text().splitlines()
    if not lines:
        return rows
    header = lines[0].split("\t")
    for line in lines[1:]:
        if line:
            rows.append(dict(zip(header, line.split("\t"))))
    return rows


def _gradient_color(stops: list[str], frac: float):
    """Interpolate a multi-stop color gradient at ``frac`` in [0, 1]."""
    frac = min(1.0, max(0.0, frac))
    if len(stops) == 1:
        return stops[0]
    pos = frac * (len(stops) - 1)
    i = min(int(pos), len(stops) - 2)
    t = pos - i
    a = np.asarray(to_rgba(stops[i]))
    b = np.asarray(to_rgba(stops[i + 1]))
    return tuple(a + (b - a) * t)


def add_cigar_read(
    ax,
    start: int,
    cigar: list[tuple[int, int]],
    length: int,
    radius: float,
    min_indel: int,
    gradient: list[str] = REDWOOD_GRADIENT,
    alpha: float = 1.0,
) -> None:
    """Draw one regular read as a CIGAR-aware, gradient-colored arc.

    Each CIGAR operation is an arc segment whose radial width encodes the
    operation — matches at the nominal width, deletions/skips thin (a
    near-gap), insertions a fat bulge back-offset by half their length — and
    whose color follows ``gradient`` from the read's genomic start to its end,
    so read direction is visible. Match ops are subdivided so the gradient
    stays smooth. Insertions/deletions shorter than ``min_indel`` bp are drawn
    as plain match, suppressing small HiFi indel noise. The read is clipped at
    one full circle so its 3' end never overlaps its 5' end.
    """
    centerline = radius - READ_ARC_WIDTH / 2
    sink: list = []                 # all wedges of this read become one PatchCollection (same look, far fewer artists)
    total = max(1, min(sum(n for op, n in cigar if op in (0, 2, 3, 7, 8)), length))
    chunk = max(1, length // 360)   # op subdivision for a smooth color gradient
    ref = start
    drawn = 0  # reference bp drawn so far — clip the read at one full circle
    for op, oplen in cigar:
        if op == 1:  # I — insertion: consumes the read, not the reference
            if oplen >= min_indel and drawn < length:
                w = CIGAR_OP_WIDTH[1]
                s = (ref - oplen // 2) % length
                color = _gradient_color(gradient, drawn / total)
                add_arc(ax, s, s + oplen, length, centerline + w / 2, w, sink=sink,
                        color=color, alpha=alpha, linewidth=0)
            continue
        if op in (4, 5, 6):  # S / H / P — no span on the circle
            continue
        # M(0) D(2) N(3) =(7) X(8) all consume the reference
        if drawn >= length:
            break  # already drew one full circle; do not overlap 3' onto 5'
        seg = min(oplen, length - drawn)
        if op in (2, 3) and oplen >= min_indel:
            w = CIGAR_OP_WIDTH[op]
        else:
            w = READ_ARC_WIDTH
        pos = 0
        while pos < seg:
            sub = min(chunk, seg - pos)
            color = _gradient_color(gradient, (drawn + pos + sub / 2) / total)
            s = (ref + pos) % length
            # Non-final sub-arcs are drawn double-length so the opaque next
            # sub-arc fully covers the anti-aliased seam between them — without
            # this, the abutting Wedge patches leave faint radial hairlines.
            draw_len = sub if pos + sub >= seg else sub * 2
            add_arc(ax, s, s + draw_len, length, centerline + w / 2, w, sink=sink,
                    color=color, alpha=alpha, linewidth=0)
            pos += sub
        ref += oplen
        drawn += seg
    if sink:
        ax.add_collection(PatchCollection(sink, match_original=True), autolim=False)


def read_reference(path: Path) -> str:
    seq = []
    for line in path.read_text().splitlines():
        if not line or line.startswith(">"):
            continue
        seq.append(line.strip())
    return "".join(seq).upper()


def sequence_length_from_gff(path: Path) -> int | None:
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            continue
        if fields[2] == "region" and "Is_circular=true" in fields[8]:
            return int(fields[4])
    return None


def sequence_length_from_bam(path: Path) -> int:
    with pysam.AlignmentFile(path, "rb") as bam:
        if not bam.lengths:
            raise ValueError(f"no reference lengths found in {path}")
        return int(bam.lengths[0])


def at_profile(reference: str, window: int = 201) -> list[float]:
    """AT fraction in a window centered on every position; windows wrap around the origin of the circular molecule."""
    n = len(reference)
    if n == 0:
        return []
    half = window // 2
    at = np.fromiter((c in "ATat" for c in reference), dtype=float, count=n)
    padded = np.take(at, np.arange(-half, n + half), mode="wrap")
    csum = np.concatenate([[0.0], np.cumsum(padded)])
    return list((csum[window:window + n] - csum[:n]) / window)


def strand_depth_profiles(bam_path: Path | None, length: int) -> tuple[np.ndarray, np.ndarray]:
    forward = np.zeros(length, dtype=float)
    reverse = np.zeros(length, dtype=float)
    if bam_path is None or not bam_path.exists():
        return forward, reverse
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            depth = reverse if read.is_reverse else forward
            for start, stop in read.get_blocks():
                span = stop - start
                if span >= length:
                    depth += 1
                    continue
                start %= length
                stop = start + span
                if stop <= length:
                    depth[start:stop] += 1
                else:
                    depth[start:length] += 1
                    depth[0 : stop - length] += 1
    return forward, reverse


def add_at_track(ax, reference: str, inner_radius: float, outer_radius: float, image_size: int = 2400) -> None:
    length = len(reference)
    if length == 0:
        return
    at_values = np.asarray(at_profile(reference))
    # Dynamic range: stretch the AT% spread (robust 2nd-98th percentiles) into
    # AT_RANGE, a sub-range of the colormap. The stretch keeps subtle
    # differences visible; the sub-range keeps the lightest values off
    # near-white and the contrast off the harsh palette extremes.
    lo, hi = np.percentile(at_values, [2, 98])
    if hi > lo:
        norm = np.clip((at_values - lo) / (hi - lo), 0.0, 1.0)
        at_values = AT_RANGE[0] + norm * (AT_RANGE[1] - AT_RANGE[0])
    else:
        at_values = np.full_like(at_values, sum(AT_RANGE) / 2)
    grid = np.zeros((image_size, image_size, 4), dtype=float)
    axis = np.linspace(-PLOT_LIMIT, PLOT_LIMIT, image_size)
    x, y = np.meshgrid(axis, axis)
    radii = np.sqrt((x * x) + (y * y))
    angles = (np.degrees(np.arctan2(y, x)) + 360) % 360
    positions = (((90 - angles) % 360) / 360 * length).astype(int) % length
    at_mask = (radii >= inner_radius) & (radii < outer_radius)
    at_colors = AT_COLORMAP(np.take(at_values, positions))
    at_colors[..., 3] = 0.92
    grid[at_mask] = at_colors[at_mask]
    ax.imshow(
        grid,
        extent=(-PLOT_LIMIT, PLOT_LIMIT, -PLOT_LIMIT, PLOT_LIMIT),
        origin="lower",
        interpolation="nearest",
        zorder=1,
    )


def add_rnaseq_depth_track(
    ax,
    forward_depth: np.ndarray,
    reverse_depth: np.ndarray,
    length: int,
    start_radius: float,
    track_width: float,
    forward_color: str,
    reverse_color: str,
    style: str = "coverage",
) -> None:
    if len(forward_depth) == 0 or len(reverse_depth) == 0:
        return
    total_depth = forward_depth + reverse_depth
    if np.max(total_depth) <= 0:
        return
    transformed_total = np.log1p(total_depth)
    max_value = float(np.max(transformed_total))
    image_size = 900
    grid = np.zeros((image_size, image_size, 4), dtype=float)
    axis = np.linspace(-PLOT_LIMIT, PLOT_LIMIT, image_size)
    x, y = np.meshgrid(axis, axis)
    radii = np.sqrt((x * x) + (y * y))
    angles = (np.degrees(np.arctan2(y, x)) + 360) % 360
    positions = (((90 - angles) % 360) / 360 * length).astype(int) % length
    heights = track_width * (np.take(transformed_total, positions) / max_value)
    band_mask = (radii >= start_radius) & (radii <= (start_radius + heights))
    if style == "strand":
        forward = np.take(forward_depth, positions)
        reverse = np.take(reverse_depth, positions)
        totals = forward + reverse
        fraction = np.divide(reverse, totals, out=np.zeros_like(reverse), where=totals > 0)
        forward_rgba = np.asarray(to_rgba(forward_color, alpha=0.82))
        reverse_rgba = np.asarray(to_rgba(reverse_color, alpha=0.82))
        colors = (forward_rgba * (1 - fraction[..., None])) + (reverse_rgba * fraction[..., None])
        grid[band_mask] = colors[band_mask]
    else:
        grid[band_mask] = np.asarray(to_rgba(forward_color, alpha=0.82))
    ax.imshow(
        grid,
        extent=(-PLOT_LIMIT, PLOT_LIMIT, -PLOT_LIMIT, PLOT_LIMIT),
        origin="lower",
        interpolation="nearest",
        zorder=2,
    )


_NAME_SUFFIX = re.compile(r"\s+(CDS|gene|mRNA|exon|rRNA|tRNA)\b.*$", re.IGNORECASE)
_NAME_PARTIAL = re.compile(r"\s*\(?\s*(?:[35]'\s*)?partial\b.*$", re.IGNORECASE)
_NAME_PAREN = re.compile(r"\s*\([^()]*\)\s*$")
_AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Sec": "U", "Pyl": "O",
}


def clean_feature_name(name: str) -> str:
    """Strip annotator decorations (MitoFinder writes 'COX1 CDS 3\' Partial CDS', 'tRNA-Asn gene', ...) and a trailing
    parenthesized comment, so a descriptive Name never becomes a label that runs off the page."""
    name = name.strip()
    name = _NAME_PAREN.sub("", name)          # "ND4L-2 gene (second copy of the duplication)" -> "ND4L-2 gene"
    name = _NAME_PARTIAL.sub("", name)
    name = _NAME_SUFFIX.sub("", name)
    return name.strip() or "feature"


def trna_short_label(name: str) -> str:
    """'tRNA-Ser2' -> 'S2', 'tRNA-Leu(UUR)' -> 'L', 'trnK' -> 'K', 'tRNA_Asn' -> 'N'; unknown names pass through."""
    match = re.search(r"(?:tRNA|trn)[-_ ]?([A-Z][a-z]{2})(\d*)", name)
    if match and match.group(1) in _AA3:
        return _AA3[match.group(1)] + match.group(2)
    match = re.search(r"^trn([A-Z])(\d*)", name)
    if match:
        return match.group(1) + match.group(2)
    return name


def collapse_locus_features(features: list[dict[str, object]]) -> list[dict[str, object]]:
    """GFF3 describes one locus with several records (gene + mRNA + CDS, gene + tRNA, ...).
    Keep the most specific record per locus (CDS / tRNA / rRNA), and drop a 'gene'/'mRNA'/'exon'
    record when a specific record with the same name overlaps it or shares its coordinates."""
    specific = [f for f in features if str(f["type"]) in {"CDS", "tRNA", "rRNA"}]
    kept = list(specific)
    for feature in features:
        ftype = str(feature["type"])
        if ftype in {"CDS", "tRNA", "rRNA"}:
            continue
        if ftype == "exon":
            continue
        start, stop, name = int(feature["start"]), int(feature["stop"]), str(feature["name"]).lower()
        covered = any(
            (int(g["start"]) == start and int(g["stop"]) == stop)
            or (str(g["name"]).lower() == name and int(g["start"]) < stop and start < int(g["stop"]))
            for g in specific
        )
        if not covered:
            kept.append(feature)
    return kept


def parse_gff(path: Path | None) -> list[dict[str, object]]:
    if path is None or not path.exists():
        return []
    features = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            continue
        _, _, feat_type, start, stop, _, strand, _, attrs = fields
        if feat_type in {"region", "source"}:
            continue
        name = feat_type
        for key in ("Name=", "gene=", "product=", "ID="):
            for attr in attrs.split(";"):
                if attr.startswith(key):
                    name = attr.split("=", 1)[1]
                    break
            else:
                continue
            break
        features.append(
            {
                "type": feat_type,
                "start": int(start) - 1,
                "stop": int(stop),
                "strand": strand,
                "name": clean_feature_name(name),
            }
        )
    return collapse_locus_features(features)


def assign_annotation_lanes(features: list[dict[str, object]]) -> list[dict[str, object]]:
    annotated = []
    lane_stops: list[int] = []
    for feature in sorted(features, key=lambda item: (int(item["start"]), -int(item["stop"]))):
        item = dict(feature)
        if item["type"] == "tRNA":
            item["lane"] = 1
            annotated.append(item)
            continue
        start = int(item["start"])
        stop = int(item["stop"])
        lane = 0
        while lane < len(lane_stops) and start < lane_stops[lane]:
            lane += 1
        if lane == len(lane_stops):
            lane_stops.append(stop)
        else:
            lane_stops[lane] = stop
        item["lane"] = min(lane, 1)
        annotated.append(item)
    return annotated


def units_per_point(ax) -> float:
    """Data units per typographic point along the x axis (equal aspect, so the same for y)."""
    fig = ax.figure
    width_in = fig.get_size_inches()[0] * ax.get_position().width
    x0, x1 = ax.get_xlim()
    return abs(x1 - x0) / (width_in * 72.0)


def label_width_units(ax, text: str, fontsize: float) -> float:
    return len(text) * fontsize * 0.66 * units_per_point(ax)


def add_feature_label(
    ax,
    feature: dict[str, object],
    length: int,
    radius: float,
    color: str,
    *,
    outer_radius: float | None = None,
    outer_color: str | None = None,
    label_text: str | None = None,
    fontsize: float = 5.0,
    min_fontsize: float = 3.6,
    prefer_outside: bool = False,
) -> None:
    """Label a feature. The name is written along the arc when it fits (shrinking the font down to
    ``min_fontsize`` first); otherwise, when ``outer_radius`` is given, it is written radially just
    outside the outer track with a short leader in the feature color, so short genes and tRNAs
    keep their labels instead of being dropped. ``prefer_outside`` skips the along-arc attempt
    (used for tRNAs, whose arcs are too thin to carry readable text)."""
    start = int(feature["start"])
    stop = int(feature["stop"])
    span = stop - start
    name = label_text if label_text is not None else str(feature["name"])
    if not name or span <= 0:
        return
    angle = theta(int(start + (span / 2)), length)
    upright = angle % 360                  # theta() runs from 90 down to -270; compare on 0-360
    rotation = angle - 90
    if 180 < upright < 360:                # along-arc text on the bottom half would read upside down
        rotation += 180
    arc_units = (span / length) * 2 * np.pi * radius
    size = fontsize if not (prefer_outside and outer_radius is not None) else min_fontsize - 1
    while size >= min_fontsize:
        if label_width_units(ax, name, size) <= arc_units * 0.92:
            x, y = polar_xy(radius, angle)
            ax.text(x, y, name, ha="center", va="center", color=color, fontsize=size, rotation=rotation,
                    rotation_mode="anchor", fontweight="bold", zorder=5)
            return
        size -= 0.4
    if outer_radius is None:
        return
    # Radial label outside the outer track, reading outward; leader from the track edge.
    # Neighboring labels (e.g. tRNA clusters) are staggered outward so they do not overprint.
    size = max(min_fontsize, fontsize - 0.8)
    placed = getattr(ax, "_redwood_outer_labels", None)
    if placed is None:
        placed = []
        ax._redwood_outer_labels = placed
    upp = units_per_point(ax)
    min_sep = size * upp * 1.3 / max(outer_radius, 1e-6) * (180 / np.pi)  # degrees
    width = label_width_units(ax, name, size)
    step = width + 0.012
    label_angle = angle
    shift_dir = 0.0                         # set once the label has been rotated off a position label
    for _ in range(12):
        tier = 0
        while any(abs(((label_angle - a + 180) % 360) - 180) < min_sep and t == tier for a, t in placed):
            tier += 1
        if tier and shift_dir:
            # already rotated off a position label and another label sits there: spread sideways instead of stacking
            # outward along the same ray (which would read as one word, e.g. "WR")
            label_angle += shift_dir * min_sep
            continue
        base = outer_radius + tier * step
        # a position label ("5,000 bp") in the way: rotate this label just clear of it; the leader bends back to the feature
        clash = _position_label_conflict(ax, label_angle, base + 0.008, base + 0.008 + width, 0.6 * size * upp)
        if clash is None:
            break
        pa, need, d = clash
        shift_dir = shift_dir or (1.0 if d >= 0 else -1.0)
        label_angle = pa + (need + 0.3) * shift_dir
    placed.append((label_angle, tier))
    lead0 = polar_xy(outer_radius - 0.010, angle)
    lead1 = polar_xy(base + 0.004, label_angle)
    ax.plot([lead0[0], lead1[0]], [lead0[1], lead1[1]], color=str(feature.get("color", color)), lw=0.6, alpha=0.9, zorder=5)
    x, y = polar_xy(base + 0.008, label_angle)
    upright = label_angle % 360
    left_half = 90 < upright < 270         # radial text on the left half would read upside down
    radial = label_angle + 180 if left_half else label_angle
    ha = "right" if left_half else "left"
    ax.text(x, y, name, ha=ha, va="center", color=outer_color or color, fontsize=size,
            rotation=radial, rotation_mode="anchor", fontweight="bold", zorder=5)


def choose_position_label_step(length: int, max_degrees: float = 60.0) -> int:
    max_bp = max(1, int(length * (max_degrees / 360)))
    nice_steps = [
        100,
        200,
        250,
        500,
        1000,
        2000,
        2500,
        5000,
        10000,
        20000,
        25000,
        50000,
        100000,
        200000,
        250000,
        500000,
        1000000,
    ]
    candidates = [step for step in nice_steps if step <= max_bp]
    return max(candidates) if candidates else nice_steps[0]


def _position_label_conflict(ax, angle: float, r0: float, r1: float, half_height: float):
    """A recorded position label that a radial label spanning radii [r0, r1] at ``angle`` would overlap, as
    (label angle, required angular separation in degrees, signed offset), or None."""
    for pa, half_deg, pr0, pr1 in getattr(ax, "_redwood_position_labels", []):
        if r1 < pr0 or r0 > pr1:
            continue
        r_mid = max(1e-6, (max(r0, pr0) + min(r1, pr1)) / 2)
        need = half_deg + float(np.degrees(half_height / r_mid)) + 0.8
        d = ((angle - pa + 180) % 360) - 180
        if abs(d) < need:
            return pa, need, d
    return None


def add_position_labels(ax, length: int, color: str) -> None:
    tick_step = 1000 if length <= 30000 else 2500
    label_step = choose_position_label_step(length)
    for bp in range(0, length, tick_step):
        angle = theta(bp, length)
        major = bp % label_step == 0
        outer_radius = 1.205 if major else 1.192
        x0, y0 = polar_xy(1.172, angle)
        x1, y1 = polar_xy(outer_radius, angle)
        ax.plot([x0, x1], [y0, y1], color=color, lw=0.7 if major else 0.55, alpha=0.62 if major else 0.36)
        if major and bp and (length - bp) < 0.35 * label_step:
            continue  # a "N,000 bp" label right next to "0 bp" overprints it
        if major:
            x, y = polar_xy(1.238, angle)
            rotation = angle - 90
            if 180 < angle % 360 < 360:
                rotation += 180
            ax.text(
                x,
                y,
                f"{bp:,} bp",
                ha="center",
                va="center",
                color=color,
                fontsize=5.2,
                rotation=rotation,
                rotation_mode="anchor",
                alpha=0.86,
            )
            # remember the label's footprint (angular half-width at its radius, radial band) so outer feature labels avoid it
            upp = units_per_point(ax)
            text_w = len(f"{bp:,} bp") * 5.2 * 0.58 * upp
            text_h = 5.2 * upp
            footprint = (angle, float(np.degrees((text_w / 2 + 0.004) / 1.238)), 1.238 - 0.6 * text_h, 1.238 + 0.6 * text_h)
            ax.__dict__.setdefault("_redwood_position_labels", []).append(footprint)


def read_spans(path: Path | None, true_length: int, max_reads: int) -> list[tuple[int, int]]:
    if path is None or not path.exists():
        return []
    spans = []
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            start = read.reference_start % true_length
            ref_len = read.query_alignment_length or read.reference_length or read.query_length or 1
            stop = start + min(ref_len, true_length)
            spans.append((start, stop))
            if len(spans) >= max_reads:
                break
    return spans


def infer_length(reference: str | Path | None, gff: str | Path | None, main_bam: str | Path | None) -> int:
    if reference:
        return len(read_reference(Path(reference)))
    if gff:
        length = sequence_length_from_gff(Path(gff))
        if length:
            return length
    if main_bam:
        return sequence_length_from_bam(Path(main_bam))
    raise ValueError("plotting requires --mito-fasta, --gff with circular region, or --main-bam")


# Ring key geometry, in data units of the map (ring radii run to ~1.18, labels to ~1.27, axes to +/-PLOT_LIMIT).
KEY_RADIUS = 0.40          # outer radius of the 90-degree cut-out (standalone maps)
KEY_RADIUS_MAX = 0.50      # on small maps the key may grow up to this radius, so its rows keep legible spacing
KEY_PITCH_PT = 5.0         # target row pitch in points when the key grows
KEY_INNER = 0.08           # inner radius (the innermost band starts here)
KEY_CLEARANCE = 1.30       # nearest point of the cut-out to the map center (clear of rings, ticks and outer labels)
KEY_LABEL_CHARS = 17       # longest key label, used to reserve room for the label column
KEY_FONT_RANGE = (3.6, 5.0)


def _key_center(radius: float = KEY_RADIUS) -> tuple[float, float]:
    """Center of the cut-out: at the bottom edge of the plot, as far left as the clearance from the map allows."""
    cy = -PLOT_LIMIT + 0.03
    cx = float(np.sqrt(max(0.0, (KEY_CLEARANCE + radius) ** 2 - cy ** 2)))
    return cx, cy


def _key_radius(n_layers: int, upp: float) -> float:
    """Outer radius of the key: KEY_RADIUS, grown towards KEY_RADIUS_MAX when the map is small (large data units per point)
    so each row gets about KEY_PITCH_PT points."""
    return float(min(KEY_RADIUS_MAX, max(KEY_RADIUS, KEY_INNER + max(1, n_layers) * KEY_PITCH_PT * upp)))


def _units_per_point_for(ax, xspan: float, yspan: float) -> float:
    """Data units per point for an equal-aspect axes showing ``xspan`` x ``yspan`` inside its allocated box
    (the original, layout-allocated position, not the aspect-adjusted one)."""
    fig = ax.figure
    pos = ax.get_position(original=True)
    w_in = fig.get_size_inches()[0] * pos.width
    h_in = fig.get_size_inches()[1] * pos.height
    return max(xspan / (w_in * 72.0), yspan / (h_in * 72.0))


def _key_font(pitch: float, upp: float) -> float:
    return float(min(KEY_FONT_RANGE[1], max(KEY_FONT_RANGE[0], 0.8 * pitch / upp)))


def track_legend_width(ax, n_layers: int = 1) -> float:
    """Extra x range (data units, right of +PLOT_LIMIT) that the key's label column needs on this axes.

    ``n_layers`` is an upper bound on the layers the key will list: fewer layers mean a larger pitch and so a larger font,
    so the reservation is made for the font of ``n_layers`` and :func:`add_track_legend` never exceeds it (stored on the
    axes). The default of one layer reserves room for the largest font."""
    width, font, radius = 0.3, KEY_FONT_RANGE[1], KEY_RADIUS
    for _ in range(8):                                   # the scale depends on the width; converges in a few steps
        upp = _units_per_point_for(ax, 2 * PLOT_LIMIT + width, 2 * PLOT_LIMIT)
        radius = _key_radius(n_layers, upp)
        cx, _ = _key_center(radius)
        font = _key_font((radius - KEY_INNER) / max(1, n_layers), upp)
        right = cx + 0.03 + KEY_LABEL_CHARS * font * 0.56 * upp + 0.02
        width = max(0.0, right - PLOT_LIMIT)
    ax._redwood_key_font = font
    ax._redwood_key_radius = radius
    return width


def track_legend_layers(*, has_rnaseq: bool, has_at: bool, has_variants: bool, has_numts: bool, has_multipass: bool,
                        has_regular: bool, has_annotation: bool = False, has_genes: bool | None = None,
                        has_trna: bool | None = None) -> list[tuple[str, str]]:
    """Layers present on a circular plot, innermost first, as ``(key, label)``. ``has_genes`` (CDS/rRNA) and ``has_trna``
    default to ``has_annotation`` when not given."""
    has_genes = has_annotation if has_genes is None else has_genes
    has_trna = has_annotation if has_trna is None else has_trna
    layers = []
    if has_regular: layers.append(("regular", "single-pass reads"))
    if has_multipass: layers.append(("multipass", "multi-pass reads"))
    if has_numts: layers.append(("numts", "NUMT loci"))
    if has_at: layers.append(("at", "AT content"))
    if has_variants: layers.append(("variants", "variant columns"))
    if has_genes:
        layers.append(("genes", "CDS / rRNA genes"))
    if has_trna:
        layers.append(("trna", "tRNA genes"))
    if has_rnaseq: layers.append(("rnaseq", "RNA-seq depth"))
    return layers


def add_track_legend(ax, layers: list[tuple[str, str]], *, dark: bool = False, has_marks: bool = True) -> int:
    """Draw a 90-degree cut-out of the ring stack as a key in the bottom-right corner of the map.

    The cut-out opens towards the map; its straight vertical edge faces a column of labels, so every band ends next to
    its own name. Bands are schematic, in the real radial order with the real colors: wood-colored single-pass reads with
    indel width changes and mismatch dots, a multi-pass read as a spiral stepping down one rung per turn, viridis NUMT
    arcs, the AT color ramp, variant-ring bars, gene blocks with the strand arrow, tRNA blocks and the RNA-seq depth
    profile. Text is sized to the band pitch (3.6-5 pt). Returns the number of layers drawn."""
    if not layers:
        return 0
    from .multipass import MULTIPASS_COLORS
    fg = "#eef4fb" if dark else "#111827"
    tick = "#9aa8b7" if dark else "#667085"
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    upp = _units_per_point_for(ax, abs(x1 - x0), abs(y1 - y0))
    n = len(layers)
    radius = getattr(ax, "_redwood_key_radius", None) or _key_radius(n, upp)
    center = _key_center(radius)
    cx, cy = center
    pitch = (radius - KEY_INNER) / n
    band = 0.78 * pitch
    fontsize = min(_key_font(pitch, upp), getattr(ax, "_redwood_key_font", KEY_FONT_RANGE[1]))
    t1, t2 = 91.0, 179.0                      # angular extent of the cut-out (degrees, counter-clockwise from +x)
    rng = np.random.default_rng(7)
    smooth = np.convolve(rng.random(160), np.ones(9) / 9, mode="same")

    def wedge(r_out, r_in, a1, a2, **kw):
        _add_patch_fast(ax, Wedge(center, r_out, a1, a2, width=r_out - r_in, linewidth=0, zorder=5, **kw))

    def arc_xy(r, a1, a2, k=40):
        a = np.radians(np.linspace(a1, a2, k))
        return cx + r * np.cos(a), cy + r * np.sin(a)

    def line_width(spacing):                  # never wider than ~55 % of the space between neighboring lines
        return float(max(0.25, min(0.9, 0.55 * spacing / upp)))

    for i, (key, label) in enumerate(layers):
        r_in = KEY_INNER + i * pitch; r_out = r_in + band; mid = (r_in + r_out) / 2
        if key == "regular":
            spacing = band / 3; lw = line_width(spacing)
            radii = [r_out - spacing * (0.5 + j) for j in range(3)]
            for r, (a1, a2), col in zip(radii, ((93, 138), (144, 177), (100, 170)), REDWOOD_GRADIENT):
                x, y = arc_xy(r, a1, a2); ax.plot(x, y, color=col, lw=lw, solid_capstyle="butt", zorder=5)
            x, y = arc_xy(radii[0], 110, 117); ax.plot(x, y, color=REDWOOD_GRADIENT[0], lw=lw * 1.5, solid_capstyle="butt", zorder=6)  # insertion
            x, y = arc_xy(radii[0], 124, 130); ax.plot(x, y, color="white", lw=lw * 0.6, solid_capstyle="butt", zorder=6)              # deletion
            if has_marks:
                for r, a, kind in ((radii[0], 100, "T"), (radii[1], 156, "A"), (radii[2], 140, "I")):
                    x, y = arc_xy(r, a, a, k=1); ax.plot(x, y, marker="o", ms=lw * 1.2, mew=0, color=MARK_COLORS[kind], zorder=7)
        elif key == "multipass":
            # one read spiraling inward: flat along each turn, stepping down one rung over a short ramp. Clockwise travel
            # (the map's direction) runs from the left edge of the cut-out (179 deg) to the top edge (91 deg).
            spacing = band / 3; lw = line_width(spacing)
            levels = [r_out - spacing * (0.5 + j) for j in range(3)]
            ramp_hi, ramp_lo = 150.0, 124.0
            for j in range(3):
                a = np.linspace(t2 - 1, ramp_hi, 20); r = np.full_like(a, levels[j])
                if j < 2:
                    a = np.concatenate([a, np.linspace(ramp_hi, ramp_lo, 14)[1:], np.linspace(ramp_lo, t1 + 1, 20)[1:]])
                    r = np.concatenate([r, np.linspace(levels[j], levels[j + 1], 14)[1:], np.full(19, levels[j + 1])])
                ax.plot(cx + r * np.cos(np.radians(a)), cy + r * np.sin(np.radians(a)), color=MULTIPASS_COLORS[0], lw=lw,
                        solid_capstyle="round", solid_joinstyle="round", zorder=5)
        elif key == "numts":
            from .numts import _identity_color
            lane = band / 2
            for (a1, a2, ident, l) in ((94, 118, 0.99, 0), (122, 136, 0.86, 0), (140, 176, 0.93, 0), (100, 150, 0.79, 1)):
                wedge(r_out - l * lane, r_out - (l + 0.85) * lane, a1, a2, facecolor=_identity_color(ident), alpha=0.95)
        elif key == "at":
            k = 44; edges = np.linspace(t1, t2, k + 1)
            for j in range(k):
                v = AT_RANGE[0] + (AT_RANGE[1] - AT_RANGE[0]) * smooth[j]
                wedge(r_out, r_in, edges[j], edges[j + 1] + 0.4, facecolor=AT_COLORMAP(v), alpha=0.92)
        elif key == "variants":
            x, y = arc_xy(r_in, t1, t2); ax.plot(x, y, color="#c8ced8", lw=0.3, zorder=5)
            bars = ((96, "mismatch", 1.0), (106, "minor", 0.3), (116, "insertion", 0.5), (126, "minor", 0.2), (136, "deletion", 1.0),
                    (146, "minor", 0.45), (156, "insertion", 0.7), (166, "mismatch", 1.0), (175, "minor", 0.6))
            for a, kind, frac in bars:
                wedge(r_in + band * max(0.15, frac), r_in, a - 1.2, a + 1.2, facecolor=VARIANT_RING_COLORS[kind], alpha=0.95)
        elif key == "genes":
            wedge(r_out, r_in, 136, 177, facecolor=FEATURE_COLORS["CDS"], alpha=0.95)
            a = np.radians(136.0); tip = np.radians(129.0)                 # arrowhead: genes point clockwise, like the map
            ax.add_patch(Polygon([(cx + r_out * np.cos(a), cy + r_out * np.sin(a)), (cx + mid * np.cos(tip), cy + mid * np.sin(tip)),
                                  (cx + r_in * np.cos(a), cy + r_in * np.sin(a))], closed=True, facecolor=FEATURE_COLORS["CDS"],
                                 linewidth=0, zorder=5))
            wedge(r_out, r_in, 94, 122, facecolor=FEATURE_COLORS["rRNA"], alpha=0.95)
        elif key == "trna":
            for a1, a2 in ((98, 104), (119, 125), (142, 148), (164, 170)):
                wedge(r_out, r_in, a1, a2, facecolor=FEATURE_COLORS["tRNA"], alpha=0.95)
        elif key == "rnaseq":
            k = 50; edges = np.linspace(t1, t2, k + 1)
            for j in range(k):
                h = band * (0.25 + 0.75 * smooth[j + 60])
                wedge(r_in + h, r_in, edges[j], edges[j + 1] + 0.4, facecolor=BARK_COLOR if (j // 7) % 2 == 0 else BARK_COLOR_ALT, alpha=0.82)
        ax.text(cx + 0.03, cy + mid, label, ha="left", va="center", fontsize=fontsize, color=fg, zorder=6)
    ax.text(cx + 0.03, cy + radius + 0.035, "ring key", ha="left", va="bottom", fontsize=fontsize, color=tick,
            fontweight="bold", zorder=6)
    return n


def draw_circular_plot(
    ax,
    *,
    length: int,
    reference: str | None = None,
    gff: Path | None = None,
    main_bam: Path | None = None,
    rnaseq_bam: Path | None = None,
    title: str | None = None,
    subtitle: str | None = None,
    max_reads: int = 80,
    dark: bool = False,
    rnaseq_style: str = "coverage",
    max_internal_gap: int = 50,
    min_pass_fraction: float = 1.0,
    wrap_ramp: float = 0.12,
    multipass: bool = True,
    min_indel: int = 10,
    trna_labels: str = "letter",
    feature_labels: bool = True,
    read_mismatches: str = "shared",
    variant_ring: bool = True,
    variant_table: Path | None = None,
    min_minor_frac: float = 0.05,
    numt_loci: Path | None = None,
    read_classes: Path | None = None,
    track_legend: bool = True,
) -> None:
    fg = "#eef4fb" if dark else "#111827"
    rna_forward = "#b6906a" if dark else BARK_COLOR
    rna_reverse = "#caa078" if dark else BARK_COLOR_ALT
    tick_color = "#9aa8b7" if dark else "#667085"
    label_color = "#ffffff"
    rnaseq_forward, rnaseq_reverse = strand_depth_profiles(rnaseq_bam, length)

    ax.set_aspect("equal")
    if track_legend:
        # upper bound on the layers the key can list, so the reserved label room is never exceeded
        key_layers = ((1 if np.max(rnaseq_forward + rnaseq_reverse, initial=0) > 0 else 0) + (2 if gff else 0) + (1 if reference else 0)
                      + (1 if variant_ring and reference and main_bam else 0) + (1 if numt_loci else 0) + (2 if main_bam else 0))
        ax.set_xlim(-PLOT_LIMIT, PLOT_LIMIT + track_legend_width(ax, max(1, key_layers)))
    else:
        ax.set_xlim(-PLOT_LIMIT, PLOT_LIMIT)
    ax.set_ylim(-PLOT_LIMIT, PLOT_LIMIT)
    ax.set_xticks([])
    ax.set_yticks([])
    layer_flags = {"has_variants": False, "has_multipass": False, "has_regular": False}

    add_position_labels(ax, length, tick_color)
    add_rnaseq_depth_track(
        ax,
        rnaseq_forward,
        rnaseq_reverse,
        length,
        1.095,
        0.079,
        rna_forward,
        rna_reverse,
        style=rnaseq_style,
    )

    has_rnaseq = rnaseq_forward is not None and np.max(rnaseq_forward + rnaseq_reverse) > 0
    outer_radius = 1.178 if has_rnaseq else 1.118
    features = assign_annotation_lanes(parse_gff(gff))
    for feature in features:
        color = FEATURE_COLORS.get(str(feature["type"]), "#d08c35")
        feature["color"] = color
        if feature["type"] == "tRNA":
            add_arc(
                ax,
                int(feature["start"]),
                int(feature["stop"]),
                length,
                1.094,
                0.030,
                color=color,
                alpha=0.95,
                linewidth=0,
            )
            if trna_labels != "none":
                text = trna_short_label(str(feature["name"])) if trna_labels == "letter" else str(feature["name"])
                add_feature_label(ax, feature, length, 1.094, label_color, outer_radius=outer_radius,
                                  outer_color=fg, label_text=text, fontsize=4.6, min_fontsize=3.8,
                                  prefer_outside=True)
            continue
        radius = 1.030 + (float(feature.get("lane", 0)) * 0.048)
        width = 0.046
        add_directional_feature(
            ax,
            int(feature["start"]),
            int(feature["stop"]),
            length,
            radius,
            width,
            str(feature["strand"]),
            color,
        )
        if feature_labels:
            add_feature_label(ax, feature, length, radius - (width / 2), label_color,
                              outer_radius=outer_radius, outer_color=fg)

    if reference:
        add_at_track(ax, reference, 0.918, 0.954)

    # Read tracks: multi-pass (rolling-circle) reads as spirals just inside the
    # annotation/AT rings, then regular single-pass reads packed inward. Per-read
    # CIGAR indel detail is drawn on both — regular reads as arc-width changes
    # (add_cigar_read), spirals as line-width changes (add_spiral_read).
    r_top, r_min, rung_w = 0.894, 0.40, 0.0087
    numt_rows = []
    if numt_loci and Path(str(numt_loci)).exists():
        from .numts import read_loci

        numt_rows = read_loci(Path(str(numt_loci)))
        add_numt_ring(ax, numt_rows, length, 0.866, 0.908)
        r_top = 0.856
    class_of: dict[str, str] = {}
    if read_classes and Path(str(read_classes)).exists():
        for line in Path(str(read_classes)).read_text().splitlines()[1:]:
            f = line.split("\t")
            if len(f) >= 2: class_of[f[0]] = f[1]
    if main_bam is not None and Path(str(main_bam)).exists():
        # Column-level disagreement of the read population with the reference (mismatch / indel table):
        # drives the variant ring and, in "shared" mode, which per-read mismatches are marked.
        variant_rows: list[dict[str, object]] = []
        if reference and (variant_ring or read_mismatches == "shared"):
            if variant_table and Path(str(variant_table)).exists():
                variant_rows = load_variant_table(Path(str(variant_table)))
            else:
                from .variants import column_variants, flagged_rows

                rows, _ = column_variants(Path(str(main_bam)), reference, min_minor_frac=min_minor_frac)
                variant_rows = flagged_rows(rows)
        mark_filter: set[int] | None = None
        if read_mismatches == "shared":
            mark_filter = {int(r["pos"]) - 1 for r in variant_rows
                           if any(tag in str(r.get("events", "")).split(",") for tag in ("mismatch", "minor"))}
        if variant_ring and reference and variant_rows:
            add_variant_ring(ax, variant_rows, length)
            layer_flags["has_variants"] = True
        want_marks = read_mismatches != "none" and reference is not None
        multipass_reads, regular_reads = classify_circular_reads(
            Path(str(main_bam)),
            length,
            max_internal_gap=max_internal_gap,
            min_pass_fraction=min_pass_fraction if multipass else 1e9,
            reference_seq=reference if want_marks else None,
            min_indel=min_indel,
        )
        # Multi-pass spirals get at most MULTIPASS_RADIUS_FRACTION of the read
        # band. classify_circular_reads returns them sorted longest-first;
        # greedily keep reads while they fit the budget, so the highest-pass
        # reads win the space and shorter ones backfill the remainder.
        total_rungs = int((r_top - r_min) / rung_w)
        multipass_budget = int(total_rungs * MULTIPASS_RADIUS_FRACTION)
        selected = []
        used = 0
        for mp in multipass_reads:
            n_rungs = max(1, int(np.ceil(mp.passes)))
            if used + n_rungs <= multipass_budget:
                selected.append((mp, n_rungs))
                used += n_rungs
        rung = 0
        for idx, (mp, n_rungs) in enumerate(selected):
            add_spiral_read(
                ax,
                mp.start,
                mp.passes,
                length,
                r_top - rung * rung_w,
                rung_w,
                wrap_ramp,
                CLASS_COLORS[class_of[mp.name]] if class_of and class_of.get(getattr(mp, "name", "")) in CLASS_COLORS and class_of.get(mp.name) not in ("mito_only", "mito_multisegment", "mito+NUMT_homology") else MULTIPASS_COLORS[idx % len(MULTIPASS_COLORS)],
                cigar=mp.cigar,
                min_indel=min_indel,
                marks=mp.marks if want_marks else None,
                mark_filter=mark_filter,
            )
            rung += n_rungs
        if selected:
            rung += 1  # blank rung separating spirals from regular reads
            layer_flags["has_multipass"] = True
        placed, _ = pack_circular_reads(regular_reads, length, pad=length * 0.004)
        for read, lane in placed:
            radius = r_top - (rung + lane) * rung_w
            if radius <= r_min:
                continue
            cls = class_of.get(getattr(read, "name", ""), None) if class_of else None
            if read.cigar:
                add_cigar_read(ax, read.start, read.cigar, length, radius,
                               min_indel, gradient=[CLASS_COLORS[cls]] if cls in CLASS_COLORS and cls not in ("mito_only", "mito_multisegment", "mito+NUMT_homology") else REDWOOD_GRADIENT)
            else:
                add_arc(
                    ax, read.start, read.start + read.span, length, radius,
                    READ_ARC_WIDTH, color=REDWOOD_GRADIENT[1], alpha=0.85,
                    linewidth=0,
                )
            if want_marks and read.marks:
                add_read_marks(ax, read.start, read.marks, length, radius, mark_filter=mark_filter)
            layer_flags["has_regular"] = True
        upp = units_per_point(ax)

        def opening_at(y):                  # usable width of the central opening at height y
            return 1.72 * float(np.sqrt(max(r_min ** 2 - y ** 2, 0.0)))

        def centered_row(items, y, size, min_size, bold):
            gap = 2.2
            width = lambda fs: sum(len(t) * fs * (0.66 if bold else 0.56) * upp for t, _ in items) + (len(items) - 1) * gap * upp
            while size > min_size and width(size) > opening_at(y):
                size -= 0.2
            x = -width(size) / 2
            for text, color in items:
                ax.text(x, y, text, ha="left", va="center", fontsize=size, color=color, zorder=6,
                        fontweight="bold" if bold else "normal")
                x += len(text) * size * (0.66 if bold else 0.56) * upp + gap * upp

        if class_of:
            centered_row([(text, CLASS_COLORS[cls]) for cls, text in (
                ("mito+nuclear_at_NUMT_locus", "NUMT junction"), ("mito+nuclear_elsewhere", "chimera"),
                ("nuclear_only_at_NUMT_locus", "nuclear"))], -0.36, 4.6, 3.6, False)
        if want_marks:
            centered_row([(t, MARK_COLORS[k]) for k, t in (("A", "A"), ("C", "C"), ("G", "G"), ("T", "T"), ("I", "ins"), ("D", "del"))],
                         -0.22, 5.2, 3.6, True)
            caption = "read mismatches" + (" (shared)" if read_mismatches == "shared" else "")
            size = 4.6
            while size > 3.6 and len(caption) * size * 0.56 * upp > opening_at(-0.29):
                size -= 0.2
            if len(caption) * size * 0.56 * upp > opening_at(-0.29):
                caption = "read mismatches"
            ax.text(0, -0.29, caption, ha="center", va="center", fontsize=size, color=tick_color, zorder=6)

    if track_legend:
        layers = track_legend_layers(has_rnaseq=bool(has_rnaseq), has_at=bool(reference),
                                     has_genes=any(str(f["type"]) != "tRNA" for f in features),
                                     has_trna=any(str(f["type"]) == "tRNA" for f in features),
                                     has_numts=bool(numt_rows), **layer_flags)
        add_track_legend(ax, layers, dark=dark, has_marks=(read_mismatches != "none" and reference is not None))
    ax.text(0, 0.02, f"{length:,}", ha="center", va="center", color=fg, fontsize=9, fontweight="bold")
    ax.text(0, -0.085, "bp", ha="center", va="center", color=tick_color, fontsize=7)
    if title:
        ax.set_title(title, color=fg, fontsize=15, fontweight="bold", pad=10)
    if subtitle:
        ax.text(
            0.5,
            -0.06,
            subtitle,
            transform=ax.transAxes,
            ha="center",
            va="top",
            color=tick_color,
            fontsize=10,
            fontstyle="italic",
        )


def plot_file(
    *,
    output_base: str,
    fileforms: list[str],
    dpi: int,
    reference_fasta: Path | None = None,
    gff: Path | None = None,
    main_bam: Path | None = None,
    rnaseq_bam: Path | None = None,
    title: str | None = None,
    subtitle: str | None = None,
    max_reads: int = 80,
    dark: bool = False,
    rnaseq_style: str = "coverage",
    transparent: bool = False,
    no_timestamp: bool = False,
    max_internal_gap: int = 50,
    min_pass_fraction: float = 1.0,
    wrap_ramp: float = 0.12,
    multipass: bool = True,
    min_indel: int = 10,
    trna_labels: str = "letter",
    feature_labels: bool = True,
    read_mismatches: str = "shared",
    variant_ring: bool = True,
    variant_table: Path | None = None,
    min_minor_frac: float = 0.05,
    numt_loci: Path | None = None,
    read_classes: Path | None = None,
    track_legend: bool = True,
) -> None:
    bg = "#0d1117" if dark else "#ffffff"
    edge = "#303946" if dark else "#d8dee8"
    reference = read_reference(reference_fasta) if reference_fasta else None
    length = infer_length(reference_fasta, gff, main_bam)
    fig, ax = plt.subplots(figsize=(5.8, 5.8), dpi=dpi)
    fig.subplots_adjust(left=0.035, right=0.965, bottom=0.055, top=0.95)
    if track_legend:       # widen the figure by the key's label column, keeping the map itself the same size
        fig.set_size_inches(5.8 * (2 * PLOT_LIMIT + track_legend_width(ax)) / (2 * PLOT_LIMIT), 5.8)
    fig.patch.set_facecolor(bg)
    ax.set_facecolor(bg)
    for spine in ax.spines.values():
        spine.set_color(edge)
        spine.set_linewidth(1.0)
    draw_circular_plot(
        ax,
        length=length,
        reference=reference,
        gff=gff,
        main_bam=main_bam,
        rnaseq_bam=rnaseq_bam,
        title=title,
        subtitle=subtitle,
        max_reads=max_reads,
        dark=dark,
        rnaseq_style=rnaseq_style,
        max_internal_gap=max_internal_gap,
        min_pass_fraction=min_pass_fraction,
        wrap_ramp=wrap_ramp,
        multipass=multipass,
        min_indel=min_indel,
        trna_labels=trna_labels,
        feature_labels=feature_labels,
        read_mismatches=read_mismatches,
        variant_ring=variant_ring,
        variant_table=variant_table,
        min_minor_frac=min_minor_frac,
        numt_loci=numt_loci,
        read_classes=read_classes,
        track_legend=track_legend,
    )
    fig.subplots_adjust(left=0.035, right=0.965, bottom=0.055, top=0.95)
    print_images(
        base_output_name=output_base,
        image_formats=fileforms,
        no_timestamp=no_timestamp,
        dpi=dpi,
        transparent=transparent,
    )
    plt.close(fig)


def run_plot(args) -> None:
    if getattr(args, "reprocess_ont", False):
        from .ont import preprocess_plot

        args = preprocess_plot(args)
    if getattr(args, "topology", "circular") == "linear":
        from .linear import run_linear_plot

        return run_linear_plot(args)
    if getattr(args, "read_classes", None):
        raise ValueError("--read-classes currently requires --topology linear")
    output_base = args.BASENAME or "redwood"
    rnaseq_style = "strand" if "rnaseq-strand" in getattr(args, "extra_tracks", []) else "coverage"
    reference_fasta = getattr(args, "mito_fasta", None) or getattr(args, "reference_fasta", None)
    title = getattr(args, "title", None)
    subtitle = getattr(args, "subtitle", None)
    plot_file(
        output_base=output_base,
        fileforms=args.fileform,
        dpi=args.dpi,
        reference_fasta=Path(reference_fasta) if reference_fasta else None,
        gff=Path(args.gff) if args.gff else None,
        main_bam=Path(args.main_bam) if args.main_bam else None,
        rnaseq_bam=Path(args.rnaseq_bam) if args.rnaseq_bam else None,
        title=title,
        subtitle=subtitle,
        max_reads=80 if getattr(args, "max_reads", None) is None else args.max_reads,
        dark=getattr(args, "dark", False),
        rnaseq_style=rnaseq_style,
        transparent=args.transparent,
        no_timestamp=args.no_timestamp,
        max_internal_gap=getattr(args, "max_internal_gap", 50),
        min_pass_fraction=getattr(args, "min_pass_fraction", 1.0),
        wrap_ramp=getattr(args, "wrap_ramp", 0.12),
        multipass=not getattr(args, "no_multipass", False),
        min_indel=getattr(args, "min_indel", 10),
        trna_labels=getattr(args, "trna_labels", "letter"),
        feature_labels=not getattr(args, "no_feature_labels", False),
        read_mismatches=getattr(args, "read_mismatches", "shared"),
        variant_ring=not getattr(args, "no_variant_ring", False),
        variant_table=Path(args.variant_table) if getattr(args, "variant_table", None) else None,
        min_minor_frac=getattr(args, "min_minor_frac", 0.05),
        numt_loci=Path(args.numt_loci) if getattr(args, "numt_loci", None) else None,
        read_classes=Path(args.circular_read_classes) if getattr(args, "circular_read_classes", None) else None,
        track_legend=not getattr(args, "no_track_legend", False),
    )


def draw_dataset_panel(ax, dataset_dir: Path, label: str, species: str, max_reads: int, dark: bool, rnaseq_style: str) -> None:
    manifest = json.loads((dataset_dir / "manifest.json").read_text())
    reference = read_reference(dataset_dir / manifest["reference"])
    draw_circular_plot(
        ax,
        length=int(manifest["sequence_length"]),
        reference=reference,
        gff=dataset_dir / manifest["annotation"],
        main_bam=dataset_dir / manifest["bam"],
        rnaseq_bam=dataset_dir / "rnaseq.mapped.bam",
        title=label,
        subtitle=species,
        max_reads=max_reads,
        dark=dark,
        rnaseq_style=rnaseq_style,
        track_legend=False,
    )
