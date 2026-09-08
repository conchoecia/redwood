"""Redwood's circular track grammar, unrolled onto a linear coordinate axis.

Palette, composition range, CIGAR widths and read gradients are shared with
the circular renderer. Layout and sequence windows never wrap the termini.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.artist import allow_rasterization
from matplotlib.collections import PolyCollection, TriMesh
from matplotlib.colors import to_hex, to_rgba
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Polygon, Rectangle
from matplotlib.textpath import TextPath
from matplotlib.tri import Triangulation
from matplotlib.transforms import Affine2D
import numpy as np

from .linear import binned, class_colors
from .linear_evidence import pack_linear_reads
from .renderer import (
    AT_COLORMAP, AT_RANGE, BARK_COLOR, BARK_COLOR_ALT, CIGAR_OP_WIDTH, FEATURE_COLORS, READ_ARC_WIDTH,
    REDWOOD_GRADIENT, _cigar_width_profile, choose_position_label_step,
)


def linear_base_fraction(sequence, composition="AT", window=201):
    """Centered base fraction, shortening windows at each terminus (no wrapping)."""
    sequence = sequence.upper()
    bases = np.frombuffer(sequence.encode("ascii"), dtype="S1")
    at = np.r_[0, np.cumsum(np.isin(bases, np.array(list(composition), dtype="S1")))]
    half = window // 2
    positions = np.arange(len(sequence))
    starts, stops = np.maximum(0, positions - half), np.minimum(len(sequence), positions + half + 1)
    return (at[stops] - at[starts]) / (stops - starts)


def annotation_rows(features):
    """Specific feature colors and separate strands, with extra rows for overlaps."""
    rows = []
    for strand in ("+", "-", ".", "?"):
        lane_stops, lanes = [], []
        genes = [f for f in features if f["type"] in {"gene", "CDS", "rRNA", "ncRNA"} and f["strand"] == strand]
        for f in sorted(genes, key=lambda f: (f["start"], -f["stop"])):
            lane = next((i for i, stop in enumerate(lane_stops) if f["start"] >= stop), len(lane_stops))
            if lane == len(lane_stops):
                lane_stops.append(f["stop"])
                lanes.append([])
            else:
                lane_stops[lane] = f["stop"]
            lanes[lane].append(f)
        rows.extend((strand, lane) for lane in lanes)
    return rows


def add_arrow(ax, start, stop, center, height, strand, length, color, alpha=.95):
    start, stop = max(.5, start + .5), min(length + .5, stop + .5)
    if stop <= start:
        return
    head = min(max(length * .0045, 70), (stop - start) * .28, 240)
    top, bottom = center - height / 2, center + height / 2
    if strand == "+":
        points = [(start, top), (stop - head, top), (stop, center), (stop - head, bottom), (start, bottom)]
    elif strand == "-":
        points = [(start, center), (start + head, top), (stop, top), (stop, bottom), (start + head, bottom)]
    else:
        points = [(start, top), (stop, top), (stop, bottom), (start, bottom)]
    ax.add_patch(Polygon(points, facecolor=color, edgecolor="none", alpha=alpha))


def read_outline(segment, center, height, min_indel, min_width=0):
    """One read silhouette, with vertices only at CIGAR width transitions."""
    profile = _cigar_width_profile(segment.cigar, min_indel, 1.0)
    widths = np.where(profile > 1, height * CIGAR_OP_WIDTH[1] / READ_ARC_WIDTH,
                      np.where(profile < 1, height * CIGAR_OP_WIDTH[2] / READ_ARC_WIDTH, height))
    widths = np.maximum(widths, min_width)
    breaks = np.r_[0, np.flatnonzero(np.diff(widths)) + 1, len(widths)]
    upper, lower = [], []
    for start, stop in zip(breaks, breaks[1:]):
        half = widths[start] / 2
        x0, x1 = segment.start + start + .5, segment.start + stop + .5
        upper.extend([(x0, center - half), (x1, center - half)])
        lower.extend([(x0, center + half), (x1, center + half)])
    return np.array(upper + lower[::-1])


class GradientRead(TriMesh):
    """Native SVG gradients with shared interpolated rendering for PDF/PNG."""

    @allow_rasterization
    def draw(self, renderer):
        from matplotlib.backends.backend_svg import RendererSVG

        # savefig can wrap the vector backend in a MixedModeRenderer.
        svg = getattr(renderer, "_vector_renderer", renderer)
        if not isinstance(svg, RendererSVG):
            return super().draw(renderer)
        if not self.get_visible():
            return
        transform = self.get_transform() + Affine2D().scale(1, -1).translate(0, svg.height)
        points = transform.transform(self.outline)
        endpoints = transform.transform([(self.outline[:, 0].min(), 0),
                                         (self.outline[:, 0].max(), 0)])
        gradient_id = f"{self.get_gid()}_gradient"
        clip_id = f"{self.get_gid()}_bounds"
        writer = svg.writer
        svg.open_group("GradientRead", gid=self.get_gid())
        writer.start("defs")
        writer.start("linearGradient", id=gradient_id, gradientUnits="userSpaceOnUse",
                     x1=str(endpoints[0, 0]), y1=str(endpoints[0, 1]),
                     x2=str(endpoints[1, 0]), y2=str(endpoints[1, 1]))
        for offset, color in zip(np.linspace(0, 1, len(self.gradient)), self.gradient):
            writer.element("stop", offset=str(offset),
                           attrib={"stop-color": to_hex(color), "stop-opacity": str(to_rgba(color)[3])})
        writer.end("linearGradient")
        bounds = self.get_clip_box()
        if bounds is not None:
            writer.start("clipPath", id=clip_id)
            writer.element("rect", x=str(bounds.x0), y=str(svg.height - bounds.y1),
                           width=str(bounds.width), height=str(bounds.height))
            writer.end("clipPath")
        writer.end("defs")
        attrs = {"fill": f"url(#{gradient_id})", "stroke": "none"}
        if bounds is not None:
            attrs["clip-path"] = f"url(#{clip_id})"
        path = "M " + " L ".join(f"{x:.6f} {y:.6f}" for x, y in points) + " Z"
        writer.element("path", d=path, attrib=attrs)
        svg.close_group("GradientRead")
        self.stale = False


def add_gradient_read(ax, segment, center, height, min_indel, gradient, min_width=0):
    """Continuous vector shading clipped to the exact CIGAR silhouette.

    Two triangles per pair of color stops replace hundreds of solid slices.
    SVG retains vector shading. The PDF exporter composites the dense read
    layer at print resolution; PNG uses the same interpolated colors.
    The wood gradient follows increasing reference coordinates, as before;
    it does not encode alignment strand or the physical ends of the read.
    """
    outline = read_outline(segment, center, height, min_indel, min_width)
    clip = Polygon(outline, closed=True, transform=ax.transData)
    if len(gradient) == 1:
        clip.set_facecolor(gradient[0])
        clip.set_edgecolor("none")
        ax.add_patch(clip)
        return clip
    x0, y0 = outline.min(axis=0)
    x1, y1 = outline.max(axis=0)
    x = np.repeat(np.linspace(x0, x1, len(gradient)), 2)
    y = np.tile([y0, y1], len(gradient))
    triangles = []
    for i in range(0, len(x) - 2, 2):
        triangles.extend([(i, i + 1, i + 2), (i + 1, i + 3, i + 2)])
    mesh = GradientRead(Triangulation(x, y, triangles),
                        facecolors=np.repeat([to_rgba(color) for color in gradient], 2, axis=0),
                        edgecolors="none", linewidths=0)
    mesh.outline, mesh.gradient = outline, gradient
    # A companion figure can contain several axes. SVG ids must be unique
    # across the whole document, including each gradient's definition.
    index = getattr(ax.figure, "_redwood_gradient_count", 0)
    ax.figure._redwood_gradient_count = index + 1
    mesh.set_gid(f"redwood_read_{index}")
    ax.add_collection(mesh)
    mesh.set_clip_path(clip)
    return mesh


def add_composition(ax, values, top, height, length):
    lo, hi = np.percentile(values, [2, 98])
    scaled = (AT_RANGE[0] + np.clip((values - lo) / (hi - lo), 0, 1) * (AT_RANGE[1] - AT_RANGE[0])
              if hi > lo else np.full_like(values, sum(AT_RANGE) / 2))
    colors = AT_COLORMAP(scaled)
    colors[:, 3] = .92
    ax.imshow(colors[np.newaxis, :, :], extent=(.5, length + .5, top + height - .035, top + .035),
              aspect="auto", interpolation="nearest")


def add_depth_band(ax, forward, reverse, top, height, length, colors, strand=False, scale="log"):
    total = forward + reverse
    maximum = float(total.max()) if len(total) else 0
    transform = np.log1p if scale == "log" else lambda x: x
    denominator = float(transform(maximum)) or 1
    heights = transform(total) / denominator * (height - .16)
    baseline = top + height - .025
    if strand:
        ratio = np.divide(reverse, total, out=np.zeros_like(reverse, dtype=float), where=total > 0)
        rgba = np.array(to_rgba(colors[0], .82))[None, :] * (1 - ratio[:, None])
        rgba += np.array(to_rgba(colors[1], .82))[None, :] * ratio[:, None]
        polygons = [[(i + .5, baseline), (i + .5, baseline - h),
                     (i + 1.5, baseline - h), (i + 1.5, baseline)] for i, h in enumerate(heights)]
        ax.add_collection(PolyCollection(polygons, facecolors=rgba, edgecolors="none"))
    else:
        x = np.arange(length + 1) + .5
        ax.fill_between(x, baseline, baseline - np.r_[heights, heights[-1]],
                        step="post", color=colors[0], alpha=.82, linewidth=0)
    return maximum


def add_endpoint_band(ax, evidence, top, height, length, bin_size, colors, clips=False):
    profiles = evidence["profiles"]
    if clips:
        left = profiles["short_left_clips"] + profiles["long_left_clips"]
        right = profiles["short_right_clips"] + profiles["long_right_clips"]
    else:
        left, right = profiles["starts"], profiles["ends"]
    x, first, widths = binned(left, bin_size)
    _, second, _ = binned(right, bin_size)
    maximum = max(int(first.max()), int(second.max()), 1)
    denominator = np.log1p(maximum)
    center, amplitude = top + height * .55, height * .32
    for side, values, sign, color in (("left", first, -1, colors[0]), ("right", second, 1, colors[1])):
        heights = sign * np.log1p(values) / denominator * amplitude
        ax.bar(x, heights, bottom=center, width=widths, align="edge", color=color, linewidth=0)
        if clips:
            _, short, _ = binned(profiles[f"short_{side}_clips"], bin_size)
            fraction = np.divide(short, values, out=np.zeros_like(short, dtype=float), where=values > 0)
            ax.bar(x, heights * fraction, bottom=center, width=widths, align="edge",
                   color=REDWOOD_GRADIENT[-1], linewidth=0)
    ax.hlines(center, .5, length + .5, color=colors[0], linewidth=.35, alpha=.35)
    return maximum


def draw_redwood_linear(args, reference, features, selected, evidence, rna):
    length = reference.length
    dark = getattr(args, "dark", False)
    fg, bg = ("#eef4fb", "#0d1117") if dark else ("#111827", "#ffffff")
    edge, muted = ("#303946", "#9aa8b7") if dark else ("#d8dee8", "#667085")
    bark = ("#b6906a", "#caa078") if dark else (BARK_COLOR, BARK_COLOR_ALT)
    requested = getattr(args, "linear_track", None)
    tracks = list(dict.fromkeys(requested if requested is not None else ["depth"]))
    if "none" in tracks:
        if len(tracks) > 1:
            raise ValueError("--linear-track none cannot be combined with other tracks")
        tracks = []
    if args.hide_evidence:
        tracks = [name for name in tracks if name not in {"ends", "clips"}]
    if getattr(args, "read_color", "wood") == "class" and not args.read_classes:
        raise ValueError("--read-color class requires --read-classes")
    placed, row_count = pack_linear_reads(selected, length)
    genes = annotation_rows(features)
    trnas = [f for f in features if f["type"] == "tRNA"]
    terminal = [f for f in features if f["type"] in {"repeat_region", "misc_feature"}]
    bands = []
    if rna is not None:
        bands.append(("rna", .50))
    if evidence and "depth" in tracks:
        bands.append(("depth", .52))
    if terminal:
        bands.append(("termini", .25))
    if trnas:
        bands.append(("trna", .18))
    bands.extend((f"genes_{i}", .23) for i in range(len(genes)))
    bands.append(("at", .20))
    if "gc" in getattr(args, "extra_tracks", []):
        bands.append(("gc", .20))
    if row_count:
        bands.append(("reads", min(2.4, max(.4, row_count * .027)) + .10))
    for name in ("ends", "clips"):
        if evidence and name in tracks:
            bands.append((name, .64))
    total = sum(h for _, h in bands) + .48
    class_mode = getattr(args, "read_color", "wood") == "class"
    class_count = len(evidence["depth"]) if evidence else 0
    legend_height = .23 * int(np.ceil(class_count / 3)) if class_mode else 0
    bottom = .54 + legend_height + (.2 if getattr(args, "show_terminal_sequences", False) else 0)
    height = total + bottom + .46
    fig = plt.figure(figsize=(args.width, height))
    fig.patch.set_facecolor(bg)
    ax = fig.add_axes([.045, bottom / height, .91, total / height])
    ax.set_facecolor(bg)
    ax.set_xlim(-length * .075, length * 1.018)
    ax.set_ylim(total, 0)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color(edge)
        spine.set_linewidth(1)
    step = choose_position_label_step(length)
    major = [1] + list(range(step, int(length * .94), step)) + [length]
    for pos in major:
        ax.vlines(pos, .24, .29, color=muted, linewidth=.55, alpha=.8)
        ax.text(pos, .17, f"{pos:,} bp", ha="left" if pos == 1 else "right" if pos == length else "center",
                va="bottom", fontsize=6.5, color=muted)
    for pos in range(max(1, step // 5), length, max(1, step // 5)):
        if pos not in major:
            ax.vlines(pos, .25, .28, color=muted, linewidth=.4, alpha=.4)
    def label(text, top, band_height):
        ax.text(-length * .012, top + band_height / 2, text, ha="right", va="center", color=muted, fontsize=6.5)
    def note(text, top, offset=.065):
        ax.text(length, top + offset, text, ha="right", va="center", color=muted, fontsize=5.8)
    top = .36
    class_palette = class_colors(evidence["depth"] if evidence else [])
    for name, band_height in bands:
        if name in {"rna", "depth"}:
            if name == "rna":
                forward, reverse = rna
                strand = "rnaseq-strand" in getattr(args, "extra_tracks", [])
                scale = "log"  # identical encoding to the circular RNA ring
            else:
                forward = sum(evidence["depth"].values(), np.zeros(length))
                reverse, strand = np.zeros(length), False
                scale = args.depth_scale
            maximum = add_depth_band(ax, forward, reverse, top, band_height, length, bark, strand, scale)
            label(getattr(args, "rnaseq_label", "RNA depth") if name == "rna" else "Read depth", top, band_height)
            note(f"Log depth · max {maximum:,.0f}×" if scale == "log" else f"0–{maximum:,.0f}×", top)
        elif name == "termini":
            itr_caps = all((f["type"] == "repeat_region" and f["attributes"].get("rpt_type") == "inverted") or
                           (f["type"] == "misc_feature" and "cap" in f["name"].lower()) for f in terminal)
            label("ITR / cap" if itr_caps else "Features", top, band_height)
            # Repeats beneath cap glyphs; the supplied strands define the mirrored arrows.
            for f in sorted(terminal, key=lambda f: f["type"] != "repeat_region"):
                is_repeat = f["type"] == "repeat_region"
                color = BARK_COLOR_ALT if is_repeat else REDWOOD_GRADIENT[0]
                add_arrow(ax, f["start"], f["stop"], top + band_height / 2, .09 if is_repeat else .16,
                          f["strand"], length, color, .6 if is_repeat else .95)
        elif name == "trna":
            label("tRNA", top, band_height)
            for f in trnas:
                add_arrow(ax, f["start"], f["stop"], top + band_height / 2, .095,
                          f["strand"], length, FEATURE_COLORS["tRNA"])
        elif name.startswith("genes_"):
            strand, row = genes[int(name.split("_")[1])]
            label(f"Genes ({strand})", top, band_height)
            for f in row:
                start, stop = max(0, f["start"]), min(length, f["stop"])
                if stop <= start:
                    continue
                add_arrow(ax, start, stop, top + band_height / 2, .16,
                          f["strand"], length, FEATURE_COLORS.get(f["type"], "#d08c35"))
                head = min(max(length * .0045, 70), (stop - start) * .28, 240)
                label_width = TextPath((0, 0), f["name"], prop=FontProperties(size=6.5, weight="bold")).get_extents().width
                available = (stop - start - head) / length * args.width * .91 / 1.093 * 72
                center = (start + stop + (head if f["strand"] == "-" else -head)) / 2
                if len(f["name"]) <= 14 and available > label_width + 3:
                    ax.text(center, top + band_height / 2, f["name"],
                            ha="center", va="center", color="white", fontsize=6.5, fontweight="bold")
        elif name in {"at", "gc"}:
            values = linear_base_fraction(reference.sequence, "AT" if name == "at" else "GC")
            add_composition(ax, values, top, band_height, length)
            label("AT %" if name == "at" else "GC %", top, band_height)
        elif name == "reads":
            label("Reads", top, band_height)
            spacing = (band_height - .08) / max(1, row_count)
            for _, group, row in placed:
                for segment in group:
                    gradient = ([class_palette[segment.read_class]] if getattr(args, "read_color", "wood") == "class"
                                else REDWOOD_GRADIENT)
                    add_gradient_read(ax, segment, top + .04 + (row + .5) * spacing,
                                      spacing * .64, getattr(args, "min_indel", 10), gradient)
        elif name in {"ends", "clips"}:
            maximum = add_endpoint_band(ax, evidence, top, band_height, length, args.bin_size, bark, name == "clips")
            label("Starts / ends" if name == "ends" else "Soft clips", top, band_height)
            note(f"{args.bin_size} bp bins · log count · max {maximum:,}", top, offset=.10)
            detail = (f"Above: left · below: right · light: 1–{args.clip_threshold - 1} bp · dark: ≥{args.clip_threshold} bp"
                      if name == "clips" else "Above: starts · below: ends")
            ax.text(1, top + .10, detail, ha="left", va="center", color=muted, fontsize=5.8)
        top += band_height
    title = getattr(args, "title", None) or reference.name
    fig.text(.5, 1 - .24 / height, title, ha="center", va="center", color=fg, fontsize=15, fontweight="bold")
    fig.text(.5, .33 / height, f"{reference.length:,} bp", ha="center", va="center", fontsize=9, fontweight="bold", color=fg)
    subtitle = getattr(args, "subtitle", None) or reference.name
    fig.text(.5, .13 / height, subtitle, ha="center", va="center", fontsize=8, fontstyle="italic", color=muted)
    if getattr(args, "show_terminal_sequences", False):
        sequence = reference.sequence.upper()
        fig.text(.11, .55 / height, f"5′-{sequence[:30]}…", fontsize=6, family="monospace", color=muted)
        fig.text(.94, .55 / height, f"…{sequence[-30:]}-3′", ha="right", fontsize=6, family="monospace", color=muted)
    if class_mode and evidence:
        handles = [Rectangle((0, 0), 1, 1, color=class_palette[name],
                             label=f"{name} ({info['reads']:,} reads)")
                   for name, info in evidence["summary"]["classes"].items()]
        legend_y = .49 + (.2 if getattr(args, "show_terminal_sequences", False) else 0)
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, legend_y / height),
                   fontsize=6, ncol=3, frameon=False, labelcolor=muted)
    # Save the actual production configuration with the numerical evidence.
    if evidence:
        evidence["summary"]["production_style"] = {
            "tracks": tracks, "read_color": getattr(args, "read_color", "wood"),
            "read_rows": row_count, "at_window": 201, "composition_windows": "truncated at termini",
            "rnaseq_label": getattr(args, "rnaseq_label", "RNA depth") if rna is not None else None,
        }
    return fig
