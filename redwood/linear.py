"""Horizontal genome figures with independent termini and auditable evidence."""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon, Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np
import pysam

from .functions import timestamp
from .linear_evidence import (
    collect_evidence, load_annotations, load_classes, load_reference,
    pack_linear_reads, read_segments, representative_segment, select_display_segments,
    terminal_group, write_evidence,
)


FEATURE_COLORS = {"gene": "#277b55", "CDS": "#277b55", "rRNA": "#ad5840",
                  "tRNA": "#9461a9", "ncRNA": "#9461a9", "repeat_region": "#447eac",
                  "misc_feature": "#86724a"}
PALETTE = ["#277b55", "#cf7b32", "#7a61a3", "#447eac", "#ab5369", "#89853e", "#637e88"]


def class_colors(labels):
    colors = {}
    for i, label in enumerate(sorted(labels)):
        lower = label.lower()
        if lower in {"unclassified", "unknown", "uninformative"}:
            colors[label] = "#999999"
        elif lower in {"all", "mito", "mito_only", "mito_like"}:
            colors[label] = PALETTE[0]
        else:
            colors[label] = PALETTE[(i % (len(PALETTE) - 1)) + 1]
    return colors


def feature_polygon(start, stop, y, strand, length):
    """Clip overhanging features to the visible molecule without wrapping them."""
    start, stop = max(.5, start + .5), min(length + .5, stop + .5)
    height = .27
    head = min(length * .006, (stop - start) * .35)
    if strand == "-":
        return [(start, y), (start + head, y + height), (stop, y + height),
                (stop, y - height), (start + head, y - height)]
    if strand == "+":
        return [(start, y + height), (stop - head, y + height), (stop, y),
                (stop - head, y - height), (start, y - height)]
    return [(start, y - height), (start, y + height), (stop, y + height), (stop, y - height)]


def draw_annotations(ax, features, length, fg, width):
    ax.hlines(0, .5, length + .5, color=fg, linewidth=1.2)
    ax.vlines([.5, length + .5], -.2, .2, color=fg, linewidth=1.2)
    genes = [f for f in features if f["type"] in {"gene", "CDS", "tRNA", "rRNA", "ncRNA"}]
    others = [f for f in features if f["type"] in {"repeat_region", "misc_feature"}]
    extremes = [0]
    for strand in ("+", "-", ".", "?"):
        lane_ends = []
        for feature in sorted((f for f in genes if f["strand"] == strand), key=lambda f: f["start"]):
            start, stop = max(0, feature["start"]), min(length, feature["stop"])
            if stop <= start:
                continue
            name = feature["name"]
            # Reserve label width as well as the genomic span so short adjacent
            # genes (e.g. ATP8 and tRNAs) keep readable labels at paper size.
            label_half = len(name) * 3.6 * length / (width * 72 * .88) / 2
            center = min(length - label_half, max(label_half, (start + stop) / 2))
            left, right = min(start, center - label_half), max(stop, center + label_half)
            lane = 0
            while lane < len(lane_ends) and left < lane_ends[lane] + length * .006:
                lane += 1
            if lane == len(lane_ends):
                lane_ends.append(right)
            else:
                lane_ends[lane] = right
            direction = -1 if strand == "-" else 1
            y = direction * (.8 + lane * 1.15)
            ax.add_patch(Polygon(feature_polygon(start, stop, y, strand, length),
                                 facecolor=FEATURE_COLORS[feature["type"]], edgecolor="none"))
            ax.text(center, y + direction * .4, name, ha="center",
                    va="bottom" if direction > 0 else "top", fontsize=7, color=fg)
            extremes.append(y + direction * 1.05)
    repeat_top = max(extremes) + .6
    for kind in ("repeat_region", "misc_feature"):
        lane = repeat_top if kind == "repeat_region" else repeat_top + 1.0
        found = False
        for f in others:
            if f["type"] != kind:
                continue
            start, stop = max(0, f["start"]), min(length, f["stop"])
            if stop <= start:
                continue
            found = True
            ax.add_patch(Polygon(feature_polygon(start, stop, lane, f["strand"], length),
                                 facecolor=FEATURE_COLORS[kind], edgecolor="none"))
            ha = "left" if start == 0 else "right" if stop == length else "center"
            x = start + 1 if ha == "left" else stop if ha == "right" else (start + stop) / 2
            ax.text(x, lane + .37, f["name"], ha=ha, va="bottom", fontsize=7, color=fg)
        if found:
            extremes.append(lane + 1.1)
    ax.set_ylim(min(extremes) - .2, max(extremes) + .2)
    ax.set_yticks([])
    ax.set_ylabel("Annotation\n+ above / − below", fontsize=9)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_reads(ax, segments, colors, length, min_indel, fg):
    groups, row_count = pack_linear_reads(segments, length)
    lines, line_colors, deletions, insertions, clips = [], [], [], [], []
    for _, group, lane in groups:
        for s in group:
            color = colors[s.read_class]
            pos = s.start + .5
            for op, n in s.cigar:
                if op in {0, 7, 8}:
                    lines.append([(pos, lane), (pos + n, lane)])
                    line_colors.append(color)
                    pos += n
                elif op in {2, 3}:
                    if n >= min_indel or op == 3:
                        deletions.append([(pos, lane), (pos + n, lane)])
                    else:
                        lines.append([(pos, lane), (pos + n, lane)])
                        line_colors.append(color)
                    pos += n
                elif op == 1 and n >= min_indel:
                    insertions.append((pos, lane))
            if s.left_clip:
                clips.append((s.start + .5, lane))
            if s.right_clip:
                clips.append((s.stop + .5, lane))
    ax.add_collection(LineCollection(lines, colors=line_colors, linewidths=.75))
    ax.add_collection(LineCollection(deletions, colors=fg, linewidths=.5, linestyles="dotted"))
    if insertions:
        ax.scatter(*zip(*insertions), marker="|", s=8, linewidths=.5, color=fg)
    if clips:
        ax.scatter(*zip(*clips), marker=".", s=3, color=fg, zorder=4)
    ax.set_ylim(max(row_count, 1), -1)
    ax.set_yticks([])
    ax.set_ylabel(f"Reads\n{len(groups):,} shown", fontsize=9)
    ax.set_title("Soft clips: endpoint dots   ·   insertions: ticks   ·   deletions/skips: dotted gaps",
                 fontsize=7, loc="left", color=fg, pad=3)


def set_depth_scale(ax, scale):
    if scale == "log":
        ax.set_yscale("function", functions=(np.log1p, np.expm1))
        high = max(1, ax.get_ylim()[1])
        ax.set_yticks([0] + [10 ** power for power in range(int(np.log10(high)) + 1)])
    else:
        ax.yaxis.set_major_locator(MaxNLocator(4, integer=True))
    ax.set_ylim(bottom=0)


def binned(values, size):
    starts = np.arange(0, len(values), size)
    return starts + .5, np.add.reduceat(values, starts), np.minimum(size, len(values) - starts)


def endpoint_scale(ax):
    ax.set_yscale("symlog", linthresh=1)
    low, high = ax.get_ylim()
    maximum = max(abs(low), abs(high), 1)
    outer = 10 ** int(np.log10(maximum))
    inner = 10 if outer >= 100 else 1
    ticks = sorted(value for value in {-outer, -inner, 0, inner, outer} if low <= value <= high)
    ax.set_yticks(ticks)
    ax.set_ylim(low, high)


def rna_depth(path, reference):
    """Primary RNA alignment base depth, without circular wrapping or pileup caps."""
    forward, reverse = np.zeros(reference.length + 1), np.zeros(reference.length + 1)
    with pysam.AlignmentFile(path, "rb") as bam:
        if reference.name not in bam.references or bam.get_reference_length(reference.name) != reference.length:
            raise ValueError("RNA BAM reference name/length does not match the linear FASTA")
        for read in bam.fetch(reference.name):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            diff = reverse if read.is_reverse else forward
            for start, stop in read.get_blocks():
                if start < 0 or stop > reference.length:
                    raise ValueError("RNA alignment extends beyond the linear reference")
                diff[start] += 1
                diff[stop] -= 1
    return np.cumsum(forward)[:-1], np.cumsum(reverse)[:-1]


def run_linear_plot(args):
    if getattr(args, "doubled", []):
        raise ValueError("--topology linear cannot be combined with --doubled; use an undoubled BAM")
    if not getattr(args, "mito_fasta", None):
        raise ValueError("--topology linear requires --mito-fasta for reference identity and terminal sequences")
    if args.max_reads < 0 or args.width <= 0 or args.bin_size < 1 or getattr(args, "min_indel", 10) < 1:
        raise ValueError("width/bin-size/min-indel must be positive; max-reads must be nonnegative")
    if min(args.terminal_window, args.junction_window, args.clip_threshold) < 1:
        raise ValueError("terminal/junction windows and clip threshold must be positive")
    reference = load_reference(args.mito_fasta)
    features = load_annotations(args.gff, reference)
    classes = load_classes(args.read_classes)
    if args.read_classes and not args.main_bam:
        raise ValueError("--read-classes requires --main-bam")
    evidence, segments, selected = None, [], []
    if args.main_bam:
        segments, counts = read_segments(args.main_bam, reference, classes if args.read_classes else None)
        evidence = collect_evidence(segments, reference, counts, args.terminal_window,
                                    args.junction_window, args.clip_threshold)
        selected = select_display_segments(segments, reference.length, getattr(args, "query", None),
                                           getattr(args, "sort", "ALNLEN"), args.max_reads,
                                           getattr(args, "linear_read_selection", "terminal-balanced"),
                                           args.terminal_window)
        groups = defaultdict(list)
        for segment in selected:
            groups[segment.name].append(segment)
        evidence["summary"]["read_selection"] = {
            "method": getattr(args, "linear_read_selection", "terminal-balanced"),
            "ranking": "ALNLEN" if getattr(args, "sort", "ALNLEN") == "POS" else getattr(args, "sort", "ALNLEN"),
            "order": "reference start, then end; pack nonoverlapping reads",
            "terminal_window": min(args.terminal_window, max(1, reference.length // 2)),
            "terminal_groups": dict(Counter(terminal_group(representative_segment(group), reference.length,
                                                           args.terminal_window) for group in groups.values())),
            "read_ids": list(groups),
        }
    rna = rna_depth(args.rnaseq_bam, reference) if args.rnaseq_bam else None
    if getattr(args, "linear_style", "redwood") == "redwood":
        from .linear_redwood import draw_redwood_linear

        fig = draw_redwood_linear(args, reference, features, selected, evidence, rna)
        return save_linear_figure(fig, args, evidence, selected)
    panels = [("annotation", 2.0)]
    if selected:
        panels.append(("reads", 2.2))
    if evidence:
        panels.append(("depth", 1.5))
        if not args.hide_evidence:
            panels.extend([("ends", 1.1), ("clips", 1.1)])
    if rna is not None:
        panels.append(("rna", 1.2))
    bg, fg = ("#101820", "#eff3f5") if getattr(args, "dark", False) else ("white", "#253342")
    height = sum(panel_height for _, panel_height in panels) + 1.4
    fig, axes = plt.subplots(len(panels), 1, sharex=True, squeeze=False,
                             figsize=(args.width, height),
                             gridspec_kw={"height_ratios": [height for _, height in panels]})
    fig.patch.set_facecolor(bg)
    axes = dict(zip((name for name, _ in panels), axes[:, 0]))
    x = np.arange(reference.length) + 1
    colors = class_colors(evidence["depth"] if evidence else [])
    for ax in axes.values():
        ax.set_facecolor(bg)
        ax.tick_params(colors=fg, labelsize=8)
        ax.yaxis.label.set_color(fg)
        ax.spines[["top", "right"]].set_visible(False)
        for spine in ax.spines.values():
            spine.set_color(fg)
            spine.set_linewidth(.5)
        ax.grid(axis="y", alpha=.12)
        ax.set_xlim(.5, reference.length + .5)
    draw_annotations(axes["annotation"], features, reference.length, fg, args.width)
    if selected:
        draw_reads(axes["reads"], selected, colors, reference.length, getattr(args, "min_indel", 10), fg)
    if evidence:
        ax = axes["depth"]
        for label, depth in evidence["depth"].items():
            n = evidence["summary"]["classes"][label]["reads"]
            ax.plot(x, depth, color=colors[label], linewidth=.85, label=f"{label} ({n:,} reads)")
        ax.set_ylabel("Read depth", fontsize=9)
        set_depth_scale(ax, args.depth_scale)
        if colors:
            ax.legend(loc="upper center", bbox_to_anchor=(.5, 1.3), frameon=False,
                      fontsize=7, ncol=min(3, len(colors)), labelcolor=fg)
        profiles = evidence["profiles"]
        if "ends" in axes:
            ax = axes["ends"]
            for key, sign, color in (("starts", 1, "#277b55"), ("ends", -1, "#ad5840")):
                bx, values, widths = binned(profiles[key], args.bin_size)
                ax.bar(bx, values * sign, width=widths, align="edge", color=color, linewidth=0)
            endpoint_scale(ax)
            ax.set_ylabel(f"Starts + / ends −\nper {args.bin_size} bp", fontsize=9)
            ax.axhline(0, color=fg, linewidth=.4)
        if "clips" in axes:
            ax = axes["clips"]
            for side, sign in (("left", 1), ("right", -1)):
                bottom = np.zeros(reference.length)
                for category, color in (("short", "#91b6c9"), ("long", "#46657f")):
                    bx, values, widths = binned(profiles[f"{category}_{side}_clips"], args.bin_size)
                    _, base, _ = binned(bottom, args.bin_size)
                    ax.bar(bx, sign * values, bottom=sign * base, width=widths, align="edge", color=color, linewidth=0)
                    bottom += profiles[f"{category}_{side}_clips"]
            endpoint_scale(ax)
            ax.set_ylabel("Soft clips\nleft + / right −", fontsize=9)
            ax.legend(handles=[Rectangle((0, 0), 1, 1, color=c, label=t) for c, t in
                      (("#91b6c9", f"1–{args.clip_threshold - 1} bp"), ("#46657f", f"≥{args.clip_threshold} bp"))],
                      frameon=False, fontsize=7, ncol=2, loc="upper center", labelcolor=fg)
    if rna is not None:
        ax = axes["rna"]
        if "rnaseq-strand" in getattr(args, "extra_tracks", []):
            ax.plot(x, rna[0], color="#447eac", lw=.8, label="RNA +")
            ax.plot(x, rna[1], color="#ad5840", lw=.8, label="RNA −")
            ax.legend(frameon=False, fontsize=7, labelcolor=fg)
        else:
            ax.fill_between(x, rna[0] + rna[1], color="#86724a", linewidth=0)
        ax.set_ylabel("RNA depth", fontsize=9)
        set_depth_scale(ax, args.depth_scale)
    last = list(axes.values())[-1]
    ticks = list(MaxNLocator(8, integer=True).tick_values(1, reference.length))
    ticks = [1] + [t for t in ticks if reference.length * .04 < t < reference.length * .94] + [reference.length]
    last.set_xticks(ticks)
    last.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    last.set_xlabel(f"Position on {reference.name} (bp)", color=fg, fontsize=10)
    title = getattr(args, "title", None) or f"{reference.name} · linear mitogenome"
    fig.suptitle(title, x=.085, ha="left", y=1 - .15 / height, color=fg, fontsize=15, fontweight="bold")
    subtitle = getattr(args, "subtitle", None)
    fig.text(.085, 1 - .47 / height, f"{reference.name}   |   {reference.length:,} bp" + (f"   |   {subtitle}" if subtitle else ""),
             fontsize=9, color=fg)
    sequence = reference.sequence.upper()
    fig.text(.085, 1 - .72 / height, f"Left  5′-{sequence[:30]}…", fontsize=7, family="monospace", color=fg)
    fig.text(.98, 1 - .72 / height, f"…{sequence[-30:]}-3′  Right", ha="right", fontsize=7, family="monospace", color=fg)
    footer = ""
    if evidence:
        footer = "Depth/end/clip tracks use all primary + supplementary alignments; the read display is limited separately."
        if classes:
            footer += " Classes supplied by user; fresh NUMTs may be indistinguishable."
        if evidence["summary"]["alignment_counts"]["hard_clipped"]:
            footer += " Hard clips: incomplete evidence."
        elif not evidence["summary"]["alignment_counts"]["supplementary"]:
            footer += " No supplementary alignments in input."
    fig.text(.085, .12 / height, footer, fontsize=6.3, color=fg, wrap=True)
    fig.subplots_adjust(left=.085, right=.98, top=1 - 1.0 / height, bottom=.68 / height, hspace=.5)
    save_linear_figure(fig, args, evidence, selected)


def save_linear_figure(fig, args, evidence, selected):
    base = args.BASENAME or "redwood"
    if not args.no_timestamp:
        base = f"{base}_{timestamp()}"
    Path(base).parent.mkdir(parents=True, exist_ok=True)
    try:
        with plt.rc_context({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}):
            for fmt in args.fileform:
                fig.savefig(f"{base}.{fmt}", format=fmt, dpi=args.dpi, transparent=args.transparent)
    finally:
        plt.close(fig)
    if evidence:
        evidence["summary"]["display"] = dict(reads=len({s.name for s in selected}),
                                              max_reads=args.max_reads, query=getattr(args, "query", None),
                                              style=getattr(args, "linear_style", "redwood"),
                                              extra_tracks=getattr(args, "linear_track", []))
        write_evidence(base, evidence)
        for warning in evidence["summary"]["warnings"]:
            print(f"redwood: {warning}", file=sys.stderr)
    print(f"Wrote linear plot: {base}")
