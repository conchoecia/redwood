"""Expanded reference-coordinate views of annotated linear genome termini.

Annotations are supplied by the caller. These panels do not infer caps,
repeat boundaries, or the reference positions of soft-clipped sequence.
"""

from __future__ import annotations

from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

from .linear_evidence import merged_intervals
from .renderer import BARK_COLOR, BARK_COLOR_ALT, FEATURE_COLORS, REDWOOD_GRADIENT


def _window_reads(selected, start, stop, maximum=12):
    """Choose unique, longest read spans overlapping a half-open zoom window."""
    groups = defaultdict(list)
    for segment in selected:
        groups[segment.name].append(segment)
    candidates = []
    for name, group in groups.items():
        visible = [s for s in group if s.start < stop and s.stop > start]
        if not visible:
            continue
        span = sum(b - a for a, b in merged_intervals((s.start, s.stop) for s in group))
        candidates.append((name, visible, span))
    chosen = sorted(candidates, key=lambda item: (-item[2], item[0]))[:maximum]
    return sorted(chosen, key=lambda item: (min(s.start for s in item[1]),
                                           max(s.stop for s in item[1]), item[0]))


def _annotation_name(feature):
    if feature["type"] == "repeat_region" and feature.get("attributes", {}).get("rpt_type") == "inverted":
        return "ITR"
    if feature["type"] == "misc_feature" and "cap" in feature["name"].lower():
        return "Cap"
    return feature["name"]


def _profile(ax, profiles, keys, start, stop, top, height, colors, bin_size, minimum_rule=.25):
    """Signed counts with shared log1p scaling above/below a zero baseline."""
    edges = np.arange(start, stop, bin_size)
    widths = np.minimum(bin_size, stop - edges)
    values = [np.array([np.sum(profiles[key][a:a + w]) for a, w in zip(edges, widths)])
              for key in keys]
    maximum = max((int(v.max()) for v in values if len(v)), default=0)
    baseline, amplitude = top + height / 2, height * .43
    denominator = np.log1p(maximum) or 1
    for counts, direction, color in zip(values, (-1, 1), colors):
        ax.bar(edges + .5, direction * np.log1p(counts) / denominator * amplitude,
               bottom=baseline, width=widths, align="edge", color=color, linewidth=0)
    ax.hlines(baseline, start + .5, stop + .5, color=colors[0],
              linewidth=max(.35, minimum_rule), alpha=.5)
    return maximum


def draw_terminal_details(args, reference, features, selected, evidence):
    """Return a two-terminal figure using selected reads and all-input profiles.

    ``args.width`` is the final figure width in inches. Widths below 5 inches
    stack the termini; wider canvases place them side by side. Metadata uses
    explicit zero-based half-open intervals as well as one-based labels.
    """
    # Import at call time so the main exporter's GradientRead handling also
    # covers these panels, without introducing a renderer import cycle.
    from .linear_redwood import add_gradient_read
    from .linear_publication import publication_font

    length = reference.length
    if length < 1:
        raise ValueError("terminal detail requires a nonempty reference")
    width = float(args.width)
    if width <= 0:
        raise ValueError("terminal detail width must be positive")
    window = min(500, length)
    family = publication_font()
    journal = getattr(args, "publication_journal", "nature")
    minimum_rule = 1.0 if journal == "nature-communications" else .25
    windows = [("Left terminus", 0, window), ("Right terminus", length - window, length)]
    dark = getattr(args, "dark", False)
    fg, bg = ("#eef4fb", "#0d1117") if dark else ("#27303b", "white")
    secondary = "#b2bfce" if dark else "#485466"
    bark = ("#b6906a", "#caa078") if dark else (BARK_COLOR, BARK_COLOR_ALT)
    relevant = {"gene", "CDS", "rRNA", "tRNA", "ncRNA", "repeat_region", "misc_feature"}
    panels = []
    for title, start, stop in windows:
        annotation = [f for f in features if f["type"] in relevant and f["start"] < stop and f["stop"] > start]
        annotation.sort(key=lambda f: (f["type"] not in {"repeat_region", "misc_feature"},
                                        f["type"] != "repeat_region", f["start"], f["stop"]))
        reads = _window_reads(selected, start, stop)
        panels.append((title, start, stop, annotation, reads))
    annotation_height = max(1, max(len(p[3]) for p in panels)) * .185
    max_rows = max(len(p[4]) for p in panels)
    reads_height = max(.18, max_rows * (.085 if max_rows <= 6 else .048))
    profiles = evidence.get("profiles", {}) if evidence else {}
    have_profiles = all(k in profiles for k in ("starts", "ends", "short_left_clips", "long_left_clips",
                                                "short_right_clips", "long_right_clips"))
    evidence_height = .86 if have_profiles else 0
    panel_height = .37 + annotation_height + .18 + reads_height + evidence_height + .34
    stacked = width < 5
    gap, outside = .20, .13
    height = 2 * panel_height + gap + 2 * outside if stacked else panel_height + 2 * outside
    fig = plt.figure(figsize=(width, height), facecolor=bg)
    metadata = {"reference": reference.name, "layout": "stacked" if stacked else "side-by-side",
                "journal": journal, "minimum_rule_pt": minimum_rule, "font_family": family,
                "annotation_source": "supplied annotations", "selection_pool": "main plot selected reads",
                "read_selection": "longest union of reference spans, then start/end coordinate order",
                "profile_population": "all input primary and supplementary alignments",
                "profile_scale": "log(1 + count), normalized separately within each profile and terminal",
                "clipping_coordinates": "reference-facing alignment ends; soft-clipped bases are not plotted",
                "panels": []}

    def text(ax, x, y, label, *, size=7, color=fg, **kwargs):
        return ax.text(x, y, label, fontsize=size, color=color, fontfamily=family, **kwargs)

    for index, (title, start, stop, annotations, reads) in enumerate(panels):
        if stacked:
            rect = [outside / width, (outside + (1 - index) * (panel_height + gap)) / height,
                    (width - 2 * outside) / width, panel_height / height]
        else:
            column_width = (width - 2 * outside - gap) / 2
            rect = [(outside + index * (column_width + gap)) / width, outside / height,
                    column_width / width, panel_height / height]
        ax = fig.add_axes(rect, facecolor=bg)
        ax.set_xlim(start + .5, stop + .5)
        ax.set_ylim(panel_height, 0)
        ax.set_axis_off()
        span = stop - start
        text(ax, start + .5, .06, title, weight="bold", va="top")
        text(ax, stop + .5, .07, f"{start + 1:,}–{stop:,} bp", size=6, color=secondary,
             ha="right", va="top")
        tick_positions = np.unique(np.r_[start + 1, np.arange(((start + 99) // 100) * 100, stop, 100), stop])
        tick_positions = tick_positions[(tick_positions >= start + 1) & (tick_positions <= stop)]
        # Avoid nearly coincident endpoint ticks when the coordinate range is
        # not a multiple of 100 (for example 15,867–16,366).
        tick_positions = [p for p in tick_positions if p in {start + 1, stop}
                          or min(p - start - 1, stop - p) >= span * .12]
        for pos in tick_positions:
            ax.vlines(pos, .27, .305, color=secondary, linewidth=max(.5, minimum_rule))
            text(ax, pos, .23, f"{pos:,}", size=6, color=secondary, va="bottom",
                 ha="left" if pos == start + 1 else "right" if pos == stop else "center")
        top = .35
        for feature in annotations:
            a, b = max(start, feature["start"]), min(stop, feature["stop"])
            is_repeat = feature["type"] == "repeat_region"
            color = (BARK_COLOR_ALT if is_repeat else REDWOOD_GRADIENT[0]
                     if feature["type"] == "misc_feature" else FEATURE_COLORS.get(feature["type"], "#2f8f46"))
            name = _annotation_name(feature)
            direction = {"+": "→", "-": "←"}.get(feature["strand"], "")
            text(ax, start + .5, top, f"{name} {direction}".strip(), size=7, va="top")
            ax.add_patch(Rectangle((a + .5, top + .105), b - a, .05,
                                   facecolor=color, edgecolor="none", alpha=.75 if is_repeat else 1))
            top += .185
        if not annotations:
            text(ax, start + .5, top, "No supplied features in this window", size=6, color=secondary, va="top")
        top = .37 + annotation_height
        text(ax, start + .5, top, f"Selected reads (n = {len(reads)})", va="top")
        text(ax, stop + .5, top + .015, "│ clipped end", size=6, color=secondary, ha="right", va="top")
        top += .18
        spacing = reads_height / max(1, max_rows)
        read_records = []
        for row, (name, segments, _) in enumerate(reads):
            center = top + (row + .5) * spacing
            record = {"read": name, "segments": []}
            for segment in sorted(segments, key=lambda s: (s.start, s.stop, s.supplementary)):
                add_gradient_read(ax, segment, center, max(spacing * .58, minimum_rule / 72),
                                  getattr(args, "min_indel", 10), REDWOOD_GRADIENT,
                                  min_width=minimum_rule / 72)
                record["segments"].append({"start": segment.start, "stop": segment.stop,
                                           "left_clip": segment.left_clip, "right_clip": segment.right_clip})
                for position, clipped, side in ((segment.start + .5, segment.left_clip, "left"),
                                                 (segment.stop + .5, segment.right_clip, "right")):
                    if clipped and start + .5 <= position <= stop + .5:
                        ax.vlines(position, center - spacing * .43, center + spacing * .43,
                                  color=REDWOOD_GRADIENT[0] if not dark else fg,
                                  linewidth=max(.7, minimum_rule),
                                  clip_on=False)
                        label_width = (len(f"{clipped} bp") * 3.6 + 4)
                        visible_span = (min(stop, segment.stop) - max(start, segment.start)) / span * rect[2] * width * 72
                        # Reserve room for both labels when both clipped ends
                        # are visible. A marker alone is preferable to overlap.
                        needed = label_width + (len(f"{segment.right_clip} bp") * 3.6 + 4
                                                if segment.left_clip and segment.right_clip else 0)
                        if max_rows <= 6 and visible_span > needed:
                            inward = 1 if side == "left" else -1
                            text(ax, position + inward * span * .015, center, f"{clipped} bp",
                                 size=6, color=fg, va="center", ha="left" if inward > 0 else "right",
                                 bbox={"facecolor": bg, "edgecolor": "none", "pad": .3})
            read_records.append(record)
        if not reads:
            text(ax, start + .5, top + .05, "No selected reads overlap this window", size=6,
                 color=secondary, va="top")
        top += reads_height + .06
        bin_size = min(5, max(1, int(getattr(args, "bin_size", 25))))
        if have_profiles:
            population = evidence.get("summary", {}).get("reads")
            text(ax, start + .5, top, "Alignment ends · above: starts / below: ends", size=6, va="top")
            top += .12
            maximum = _profile(ax, profiles, ("starts", "ends"), start, stop, top, .20, bark, bin_size,
                               minimum_rule)
            top += .24
            text(ax, start + .5, top, "Soft-clipped ends · above: left / below: right", size=6, va="top")
            clip_profiles = {"left": profiles["short_left_clips"] + profiles["long_left_clips"],
                             "right": profiles["short_right_clips"] + profiles["long_right_clips"]}
            top += .12
            clip_maximum = _profile(ax, clip_profiles, ("left", "right"), start, stop, top, .20,
                                    bark, bin_size, minimum_rule)
            top += .24
            text(ax, start + .5, top, f"{bin_size} bp bins · log counts · maxima {maximum:,} / {clip_maximum:,}",
                 size=6, color=secondary, va="top")
            top += .13
            if population is not None:
                text(ax, start + .5, top, f"Profiles: all {population:,} input reads", size=6,
                     color=secondary, va="top")
            top += .16
        # An explicitly reference-derived sequence label; no clipped read bases
        # are placed on this coordinate axis or presented as cap sequence.
        count = min(18, length)
        seq_start, seq_stop = (0, count) if index == 0 else (length - count, length)
        sequence = reference.sequence[seq_start:seq_stop].upper()
        sequence_label = f"Ref. {seq_start + 1:,}–{seq_stop:,}: {sequence}"
        available_points = rect[2] * width * 72
        if len(sequence_label) * 3.6 < available_points:
            text(ax, start + .5, panel_height - .08, sequence_label, size=6, color=secondary, va="bottom")
        metadata["panels"].append({"label": title, "start": start, "stop": stop,
                                    "coordinate_system": "zero-based half-open",
                                    "first_base": start + 1, "last_base": stop,
                                    "read_ids": [r[0] for r in reads], "reads": read_records,
                                    "annotations": [{k: f[k] for k in ("name", "type", "start", "stop", "strand")}
                                                    for f in annotations],
                                    "profile_bin_size": bin_size if have_profiles else None})
    fig._redwood_terminal_metadata = metadata
    shown = " and ".join(str(len(panel["read_ids"])) for panel in metadata["panels"])
    caption = (f"Terminal detail for {reference.name}: reference positions 1–{window:,} and "
               f"{length - window + 1:,}–{length:,}. Features and reference sequence come from the supplied "
               f"annotation and FASTA. The panels display {shown} unique reads, respectively, chosen by "
               "longest reference span from the main plot's selected reads and arranged by alignment "
               "coordinates. Brown shading follows reference position and does not indicate read strand. "
               "Thin read segments indicate alignment indels; their display width is bounded by the "
               f"{minimum_rule:g} pt minimum rule. Vertical markers identify soft-clipped alignment ends; "
               "unaligned clipped sequence is not placed on the reference axis. ")
    if have_profiles:
        count = evidence.get("summary", {}).get("reads")
        population = f"all {count:,} input reads" if count is not None else "the full input population"
        caption += (f"Endpoint and clip profiles use {population}, including primary and supplementary "
                    f"alignments, rather than only the displayed reads. Counts use {bin_size} bp bins and "
                    "log(1 + count), normalized independently for each profile and terminal. Starts and "
                    "left soft-clipped ends are above the baseline; ends and right soft-clipped ends "
                    "are below. Left and right refer to reference coordinates.")
    fig._redwood_caption = caption
    return fig
