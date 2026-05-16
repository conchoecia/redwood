#!/usr/bin/env python
"""Build the light/dark README preview figure from committed fixtures."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/redwood-matplotlib")

import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import numpy as np
import pysam


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


def depth_profile(bam_path: Path, length: int) -> np.ndarray:
    depth = np.zeros(length, dtype=float)
    if not bam_path.exists():
        return depth
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            for start, stop in read.get_blocks():
                start %= length
                stop = min(stop, start + length)
                if stop <= length:
                    depth[start:stop] += 1
                else:
                    depth[start:length] += 1
                    depth[0 : stop - length] += 1
    return depth


def add_at_track(ax, reference: str, inner_radius: float, outer_radius: float, image_size: int = 760) -> None:
    length = len(reference)
    at_values = np.asarray(at_profile(reference))
    grid = np.zeros((image_size, image_size, 4), dtype=float)
    axis = np.linspace(-1.22, 1.22, image_size)
    x, y = np.meshgrid(axis, axis)
    radii = np.sqrt((x * x) + (y * y))
    angles = (np.degrees(np.arctan2(y, x)) + 360) % 360
    positions = (((90 - angles) % 360) / 360 * length).astype(int) % length
    at_mask = (radii >= inner_radius) & (radii < outer_radius)
    at_colors = plt.cm.cividis(np.take(at_values, positions))
    at_colors[..., 3] = 0.92
    grid[at_mask] = at_colors[at_mask]
    ax.imshow(grid, extent=(-1.22, 1.22, -1.22, 1.22), origin="lower", interpolation="nearest", zorder=1)


def add_rnaseq_depth_track(
    ax,
    depth: np.ndarray,
    length: int,
    start_radius: float,
    track_width: float,
    color: str,
    bins: int = 900,
) -> None:
    if len(depth) == 0 or np.max(depth) <= 0:
        return
    bins = min(bins, length)
    transformed = np.log1p(depth)
    max_value = float(np.max(transformed))
    for bin_index in range(bins):
        start = int(bin_index * length / bins)
        stop = int((bin_index + 1) * length / bins)
        if stop <= start:
            continue
        value = float(np.mean(transformed[start:stop]))
        if value <= 0:
            continue
        width = track_width * (value / max_value)
        add_arc(
            ax,
            start,
            stop,
            length,
            start_radius + width,
            width,
            color=color,
            alpha=0.76,
            linewidth=0,
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


def draw_panel(ax, dataset_dir: Path, label: str, species: str, max_reads: int, dark: bool) -> None:
    manifest = json.loads((dataset_dir / "manifest.json").read_text())
    length = int(manifest["sequence_length"])
    reference = read_reference(dataset_dir / manifest["reference"])
    fg = "#eef4fb" if dark else "#111827"
    read_color = "#e9eef5" if dark else "#111111"
    rna_color = "#d8e0ea" if dark else "#111111"
    tick_color = "#9aa8b7" if dark else "#667085"
    rnaseq_bam = dataset_dir / "rnaseq.mapped.bam"
    rnaseq_depth = depth_profile(rnaseq_bam, length)

    ax.set_aspect("equal")
    ax.set_xlim(-1.22, 1.22)
    ax.set_ylim(-1.22, 1.22)
    ax.set_xticks([])
    ax.set_yticks([])

    for radius, width, alpha in [(0.24, 0.012, 0.28)]:
        ax.add_patch(Wedge((0, 0), radius, 0, 360, width=width, color=tick_color, alpha=alpha))

    for bp in range(0, length, 5000):
        angle = theta(bp, length)
        x0, y0 = 1.145 * np.cos(np.deg2rad(angle)), 1.145 * np.sin(np.deg2rad(angle))
        x1, y1 = 1.19 * np.cos(np.deg2rad(angle)), 1.19 * np.sin(np.deg2rad(angle))
        ax.plot([x0, x1], [y0, y1], color=tick_color, lw=0.8, alpha=0.65)

    add_rnaseq_depth_track(ax, rnaseq_depth, length, 1.065, 0.079, rna_color)

    for feature in parse_gff(dataset_dir / manifest["annotation"]):
        color = FEATURE_COLORS.get(str(feature["type"]), "#d08c35")
        radius = 1.030 if feature["type"] == "tRNA" else 0.986
        width = 0.032 if feature["type"] == "tRNA" else 0.050
        add_arc(ax, int(feature["start"]), int(feature["stop"]), length, radius, width, color=color, alpha=0.95, linewidth=0)

    add_at_track(ax, reference, 0.918, 0.954)

    for row, (start, stop) in enumerate(read_spans(dataset_dir / manifest["bam"], length, max_reads)):
        radius = 0.894 - (row * 0.0087)
        add_arc(ax, start, stop, length, radius, 0.0058, color=read_color, alpha=0.72, linewidth=0)

    ax.text(0, 0.02, f"{length:,}", ha="center", va="center", color=fg, fontsize=9, fontweight="bold")
    ax.text(0, -0.085, "bp", ha="center", va="center", color=tick_color, fontsize=7)
    ax.set_title(label, color=fg, fontsize=15, fontweight="bold", pad=10)
    ax.text(0.5, -0.06, species, transform=ax.transAxes, ha="center", va="top", color=tick_color, fontsize=10, fontstyle="italic")


def write_grid(datasets_dir: Path, output: Path, dark: bool) -> None:
    bg = "#0d1117" if dark else "#ffffff"
    edge = "#303946" if dark else "#d8dee8"
    fig, axes = plt.subplots(2, 2, figsize=(10, 11), dpi=170)
    fig.patch.set_facecolor(bg)
    for ax, (dataset, label, species, max_reads) in zip(axes.flat, DATASETS):
        ax.set_facecolor(bg)
        for spine in ax.spines.values():
            spine.set_color(edge)
            spine.set_linewidth(1.0)
        draw_panel(ax, datasets_dir / dataset, label, species, max_reads, dark)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(left=0.035, right=0.965, bottom=0.055, top=0.95, wspace=0.10, hspace=0.30)
    fig.savefig(output, facecolor=bg)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets-dir", type=Path, default=Path("examples/datasets"))
    parser.add_argument("--outdir", type=Path, default=Path("docs/assets"))
    args = parser.parse_args(argv)

    write_grid(args.datasets_dir, args.outdir / "redwood-grid-light.png", dark=False)
    write_grid(args.datasets_dir, args.outdir / "redwood-grid-dark.png", dark=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
