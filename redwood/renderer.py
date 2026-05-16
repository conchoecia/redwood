"""Modern redwood circular plotting backend."""

from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/redwood-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Polygon, Wedge
import numpy as np
import pysam

from .functions import print_images


PLOT_LIMIT = 1.30

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
    ax.add_patch(
        Polygon(
            [polar_xy(inner, base), polar_xy(radius, base), polar_xy(middle, tip)],
            closed=True,
            facecolor=color,
            edgecolor="none",
            alpha=alpha,
        )
    )


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


def add_at_track(ax, reference: str, inner_radius: float, outer_radius: float, image_size: int = 760) -> None:
    length = len(reference)
    if length == 0:
        return
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


def add_position_labels(ax, length: int, color: str) -> None:
    step = 1000 if length <= 30000 else 2500
    for bp in range(0, length, step):
        angle = theta(bp, length)
        major = bp % 5000 == 0
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
) -> None:
    fg = "#eef4fb" if dark else "#111827"
    read_color = "#e9eef5" if dark else "#111111"
    rna_forward = "#aeb9c5" if dark else "#3b414a"
    rna_reverse = "#d79466" if dark else "#b85f36"
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

    for row, (start, stop) in enumerate(read_spans(main_bam, length, max_reads)):
        radius = 0.894 - (row * 0.0087)
        if radius <= 0.28:
            break
        add_arc(ax, start, stop, length, radius, 0.0058, color=read_color, alpha=0.72, linewidth=0)

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
        max_reads=getattr(args, "max_reads", 80),
        dark=getattr(args, "dark", False),
        rnaseq_style=rnaseq_style,
        transparent=args.transparent,
        no_timestamp=args.no_timestamp,
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
