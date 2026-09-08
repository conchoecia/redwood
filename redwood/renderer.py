"""Modern redwood circular plotting backend."""

from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/redwood-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
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
# genomic start, through cinnamon, to pale sapwood at its end. The colour
# sweep makes read direction visible.
REDWOOD_GRADIENT = ["#3a1d10", "#9c5a2c", "#d3a878"]

# Tree-anatomy theming: the RNA-seq depth ring is recoloured as bark, and the
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


def add_arc(ax, start: int, stop: int, length: int, radius: float, width: float, **kwargs) -> None:
    if stop <= start:
        return
    if stop - start >= length:
        ax.add_patch(Wedge((0, 0), radius, 0, 360, width=width, **kwargs))
        return
    if stop > length:
        add_arc(ax, start, length, length, radius, width, **kwargs)
        add_arc(ax, 0, stop - length, length, radius, width, **kwargs)
        return
    ax.add_patch(Wedge((0, 0), radius, theta(stop, length), theta(start, length), width=width, **kwargs))


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
    matches, thin for deletions, fat for insertions (centred, back-offset by
    half their length). Indels shorter than ``min_indel`` bp keep the match
    width. Insertions are overlaid last so a following match cannot clobber
    them."""
    thin = base_lw * 0.34
    fat = base_lw * 1.75
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
    base_linewidth: float = 0.9,
    points_per_turn: int = 540,
) -> None:
    """Draw a multi-pass (rolling-circle) read as an inward spiral.

    The read traverses the circle ``passes`` times. Each completed turn drops
    one rung; the radius is flat for the first ``1 - wrap_ramp`` of a turn and
    then ramps diagonally down over the final ``wrap_ramp`` fraction, so the
    step into the next rung is visible. When a CIGAR is supplied the spiral's
    line width varies along its length — fat at insertions, thin at deletions —
    so multi-pass reads carry the same indel detail as the regular reads.
    """
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


def _gradient_color(stops: list[str], frac: float):
    """Interpolate a multi-stop colour gradient at ``frac`` in [0, 1]."""
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
    """Draw one regular read as a CIGAR-aware, gradient-coloured arc.

    Each CIGAR operation is an arc segment whose radial width encodes the
    operation — matches at the nominal width, deletions/skips thin (a
    near-gap), insertions a fat bulge back-offset by half their length — and
    whose colour follows ``gradient`` from the read's genomic start to its end,
    so read direction is visible. Match ops are subdivided so the gradient
    stays smooth. Insertions/deletions shorter than ``min_indel`` bp are drawn
    as plain match, suppressing small HiFi indel noise. The read is clipped at
    one full circle so its 3' end never overlaps its 5' end.
    """
    centerline = radius - READ_ARC_WIDTH / 2
    total = max(1, min(sum(n for op, n in cigar if op in (0, 2, 3, 7, 8)), length))
    chunk = max(1, length // 360)   # op subdivision for a smooth colour gradient
    ref = start
    drawn = 0  # reference bp drawn so far — clip the read at one full circle
    for op, oplen in cigar:
        if op == 1:  # I — insertion: consumes the read, not the reference
            if oplen >= min_indel and drawn < length:
                w = CIGAR_OP_WIDTH[1]
                s = (ref - oplen // 2) % length
                color = _gradient_color(gradient, drawn / total)
                add_arc(ax, s, s + oplen, length, centerline + w / 2, w,
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
            add_arc(ax, s, s + draw_len, length, centerline + w / 2, w,
                    color=color, alpha=alpha, linewidth=0)
            pos += sub
        ref += oplen
        drawn += seg


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
    half = window // 2
    doubled = reference + reference
    profile = []
    offset = len(reference)
    for i in range(len(reference)):
        segment = doubled[offset + i - half : offset + i + half + 1]
        if not segment:
            profile.append(0.0)
        else:
            profile.append((segment.count("A") + segment.count("T")) / len(segment))
    return profile


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
        for attr in attrs.split(";"):
            if attr.startswith("Name="):
                name = attr.split("=", 1)[1]
                break
        features.append(
            {
                "type": feat_type,
                "start": int(start) - 1,
                "stop": int(stop),
                "strand": strand,
                "name": name,
            }
        )
    return features


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


def add_feature_label(ax, feature: dict[str, object], length: int, radius: float, color: str) -> None:
    start = int(feature["start"])
    stop = int(feature["stop"])
    span = stop - start
    if feature["type"] == "tRNA" or span < 340:
        return
    name = str(feature["name"])
    if not name or len(name) > 14:
        return
    angle = theta(int(start + (span / 2)), length)
    x, y = polar_xy(radius, angle)
    rotation = angle - 90
    if 90 < angle < 270:
        rotation += 180
    ax.text(
        x,
        y,
        name,
        ha="center",
        va="center",
        color=color,
        fontsize=5.0,
        rotation=rotation,
        rotation_mode="anchor",
        fontweight="bold",
        zorder=5,
    )


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
        if major:
            x, y = polar_xy(1.238, angle)
            rotation = angle - 90
            if 90 < angle < 270:
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
) -> None:
    fg = "#eef4fb" if dark else "#111827"
    rna_forward = "#b6906a" if dark else BARK_COLOR
    rna_reverse = "#caa078" if dark else BARK_COLOR_ALT
    tick_color = "#9aa8b7" if dark else "#667085"
    label_color = "#ffffff"
    rnaseq_forward, rnaseq_reverse = strand_depth_profiles(rnaseq_bam, length)

    ax.set_aspect("equal")
    ax.set_xlim(-PLOT_LIMIT, PLOT_LIMIT)
    ax.set_ylim(-PLOT_LIMIT, PLOT_LIMIT)
    ax.set_xticks([])
    ax.set_yticks([])

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

    for feature in assign_annotation_lanes(parse_gff(gff)):
        color = FEATURE_COLORS.get(str(feature["type"]), "#d08c35")
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
        add_feature_label(ax, feature, length, radius - (width / 2), label_color)

    if reference:
        add_at_track(ax, reference, 0.918, 0.954)

    # Read tracks: multi-pass (rolling-circle) reads as spirals just inside the
    # annotation/AT rings, then regular single-pass reads packed inward. Per-read
    # CIGAR indel detail is drawn on both — regular reads as arc-width changes
    # (add_cigar_read), spirals as line-width changes (add_spiral_read).
    r_top, r_min, rung_w = 0.894, 0.40, 0.0087
    if main_bam is not None and Path(str(main_bam)).exists():
        multipass_reads, regular_reads = classify_circular_reads(
            Path(str(main_bam)),
            length,
            max_internal_gap=max_internal_gap,
            min_pass_fraction=min_pass_fraction if multipass else 1e9,
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
                MULTIPASS_COLORS[idx % len(MULTIPASS_COLORS)],
                cigar=mp.cigar,
                min_indel=min_indel,
            )
            rung += n_rungs
        if selected:
            rung += 1  # blank rung separating spirals from regular reads
        placed, _ = pack_circular_reads(regular_reads, length, pad=length * 0.004)
        for read, lane in placed:
            radius = r_top - (rung + lane) * rung_w
            if radius <= r_min:
                continue
            if read.cigar:
                add_cigar_read(ax, read.start, read.cigar, length, radius,
                               min_indel)
            else:
                add_arc(
                    ax, read.start, read.start + read.span, length, radius,
                    READ_ARC_WIDTH, color=REDWOOD_GRADIENT[1], alpha=0.85,
                    linewidth=0,
                )

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
) -> None:
    bg = "#0d1117" if dark else "#ffffff"
    edge = "#303946" if dark else "#d8dee8"
    reference = read_reference(reference_fasta) if reference_fasta else None
    length = infer_length(reference_fasta, gff, main_bam)
    fig, ax = plt.subplots(figsize=(5.8, 5.8), dpi=dpi)
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
    )
