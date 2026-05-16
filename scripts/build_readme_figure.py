#!/usr/bin/env python
"""Build the light/dark README preview figure from committed fixtures."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/redwood-matplotlib")

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Polygon, Wedge
import numpy as np
import pysam


PLOT_LIMIT = 1.30

DATASETS = [
    ("human", "Human", "Homo sapiens", 24),
    ("mouse", "Mouse", "Mus musculus", 36),
    ("drosophila", "Fly", "Drosophila melanogaster", 44),
    ("sponge", "Sponge", "Ephydatia muelleri", 56),
]

FEATURE_COLORS = {
    "gene": "#2f8f46",
    "CDS": "#2f8f46",
    "rRNA": "#c23b3b",
    "tRNA": "#d870a2",
}

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
    points = [
        polar_xy(inner, base),
        polar_xy(radius, base),
        polar_xy(middle, tip),
    ]
    ax.add_patch(Polygon(points, closed=True, facecolor=color, edgecolor="none", alpha=alpha))


def read_reference(path: Path) -> str:
    seq = []
    for line in path.read_text().splitlines():
        if not line or line.startswith(">"):
            continue
        seq.append(line.strip())
    return "".join(seq).upper()


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


def strand_depth_profiles(bam_path: Path, length: int) -> tuple[np.ndarray, np.ndarray]:
    forward = np.zeros(length, dtype=float)
    reverse = np.zeros(length, dtype=float)
    if not bam_path.exists():
        return forward, reverse
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            depth = reverse if read.is_reverse else forward
            for start, stop in read.get_blocks():
                start %= length
                stop = min(stop, start + length)
                if stop <= length:
                    depth[start:stop] += 1
                else:
                    depth[start:length] += 1
                    depth[0 : stop - length] += 1
    return forward, reverse


def add_at_track(ax, reference: str, inner_radius: float, outer_radius: float, image_size: int = 760) -> None:
    length = len(reference)
    at_values = np.asarray(at_profile(reference))
    grid = np.zeros((image_size, image_size, 4), dtype=float)
    axis = np.linspace(-PLOT_LIMIT, PLOT_LIMIT, image_size)
    x, y = np.meshgrid(axis, axis)
    radii = np.sqrt((x * x) + (y * y))
    angles = (np.degrees(np.arctan2(y, x)) + 360) % 360
    positions = (((90 - angles) % 360) / 360 * length).astype(int) % length
    at_mask = (radii >= inner_radius) & (radii < outer_radius)
    at_colors = plt.cm.cividis(np.take(at_values, positions))
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
        coverage_rgba = np.asarray(to_rgba(forward_color, alpha=0.82))
        grid[band_mask] = coverage_rgba
    ax.imshow(
        grid,
        extent=(-PLOT_LIMIT, PLOT_LIMIT, -PLOT_LIMIT, PLOT_LIMIT),
        origin="lower",
        interpolation="nearest",
        zorder=2,
    )


def parse_gff(path: Path) -> list[dict[str, object]]:
    features = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            continue
        seq, source, feat_type, start, stop, score, strand, phase, attrs = fields
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
    center = start + (span / 2)
    angle = theta(int(center), length)
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
        fontsize=4.8,
        rotation=rotation,
        rotation_mode="anchor",
        fontweight="bold",
        zorder=5,
    )


def add_position_labels(ax, length: int, color: str) -> None:
    step = 5000 if length >= 10000 else 2500
    for bp in range(0, length, step):
        angle = theta(bp, length)
        x0, y0 = polar_xy(1.172, angle)
        x1, y1 = polar_xy(1.205, angle)
        ax.plot([x0, x1], [y0, y1], color=color, lw=0.7, alpha=0.58)
        label_radius = 1.238
        x, y = polar_xy(label_radius, angle)
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


def read_spans(path: Path, true_length: int, max_reads: int) -> list[tuple[int, int]]:
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


def draw_panel(
    ax,
    dataset_dir: Path,
    label: str,
    species: str,
    max_reads: int,
    dark: bool,
    rnaseq_style: str,
) -> None:
    manifest = json.loads((dataset_dir / "manifest.json").read_text())
    length = int(manifest["sequence_length"])
    reference = read_reference(dataset_dir / manifest["reference"])
    fg = "#eef4fb" if dark else "#111827"
    read_color = "#e9eef5" if dark else "#111111"
    rna_forward = "#aeb9c5" if dark else "#3b414a"
    rna_reverse = "#d79466" if dark else "#b85f36"
    tick_color = "#9aa8b7" if dark else "#667085"
    label_color = "#f4f7fb" if dark else "#0f172a"
    rnaseq_bam = dataset_dir / "rnaseq.mapped.bam"
    rnaseq_forward, rnaseq_reverse = strand_depth_profiles(rnaseq_bam, length)

    ax.set_aspect("equal")
    ax.set_xlim(-PLOT_LIMIT, PLOT_LIMIT)
    ax.set_ylim(-PLOT_LIMIT, PLOT_LIMIT)
    ax.set_xticks([])
    ax.set_yticks([])

    for radius, width, alpha in [(0.24, 0.012, 0.28)]:
        ax.add_patch(Wedge((0, 0), radius, 0, 360, width=width, color=tick_color, alpha=alpha))

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

    features = assign_annotation_lanes(parse_gff(dataset_dir / manifest["annotation"]))
    for feature in features:
        color = FEATURE_COLORS.get(str(feature["type"]), "#d08c35")
        if feature["type"] == "tRNA":
            radius = 1.086
            width = 0.024
            add_arc(
                ax,
                int(feature["start"]),
                int(feature["stop"]),
                length,
                radius,
                width,
                color=color,
                alpha=0.95,
                linewidth=0,
            )
            continue
        radius = 1.020 + (float(feature.get("lane", 0)) * 0.038)
        width = 0.035
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

    add_at_track(ax, reference, 0.918, 0.954)

    for row, (start, stop) in enumerate(read_spans(dataset_dir / manifest["bam"], length, max_reads)):
        radius = 0.894 - (row * 0.0087)
        add_arc(ax, start, stop, length, radius, 0.0058, color=read_color, alpha=0.72, linewidth=0)

    ax.text(0, 0.02, f"{length:,}", ha="center", va="center", color=fg, fontsize=9, fontweight="bold")
    ax.text(0, -0.085, "bp", ha="center", va="center", color=tick_color, fontsize=7)
    ax.set_title(label, color=fg, fontsize=15, fontweight="bold", pad=10)
    ax.text(
        0.5,
        -0.06,
        species,
        transform=ax.transAxes,
        ha="center",
        va="top",
        color=tick_color,
        fontsize=10,
        fontstyle="italic",
    )


def write_grid(datasets_dir: Path, output: Path, dark: bool, rnaseq_style: str) -> None:
    bg = "#0d1117" if dark else "#ffffff"
    edge = "#303946" if dark else "#d8dee8"
    fig, axes = plt.subplots(2, 2, figsize=(10, 11), dpi=170)
    fig.patch.set_facecolor(bg)
    for ax, (dataset, label, species, max_reads) in zip(axes.flat, DATASETS):
        ax.set_facecolor(bg)
        for spine in ax.spines.values():
            spine.set_color(edge)
            spine.set_linewidth(1.0)
        draw_panel(ax, datasets_dir / dataset, label, species, max_reads, dark, rnaseq_style)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(left=0.035, right=0.965, bottom=0.055, top=0.95, wspace=0.10, hspace=0.30)
    fig.savefig(output, facecolor=bg)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets-dir", type=Path, default=Path("examples/datasets"))
    parser.add_argument("--outdir", type=Path, default=Path("docs/assets"))
    parser.add_argument(
        "--rnaseq-style",
        choices=["coverage", "strand"],
        default="coverage",
        help="Render RNA-seq as a plain coverage histogram or color it by alignment strand.",
    )
    args = parser.parse_args(argv)

    write_grid(args.datasets_dir, args.outdir / "redwood-grid-light.png", dark=False, rnaseq_style=args.rnaseq_style)
    write_grid(args.datasets_dir, args.outdir / "redwood-grid-dark.png", dark=True, rnaseq_style=args.rnaseq_style)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
