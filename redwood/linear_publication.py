"""Physical-size publication layouts for Redwood's linear tracks.

All vertical coordinates are points. Text stays the same size when column
width changes; feature names get external lanes instead of being discarded.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_rgba
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.patches import Rectangle
from matplotlib.textpath import TextPath
import numpy as np

from .linear import binned, class_colors
from .linear_evidence import pack_linear_reads
from .linear_redwood import add_arrow, add_gradient_read, annotation_rows, linear_base_fraction
from .renderer import AT_COLORMAP, AT_RANGE, BARK_COLOR, BARK_COLOR_ALT, FEATURE_COLORS, REDWOOD_GRADIENT


def publication_font():
    for family in ("Arial", "Helvetica"):
        try:
            findfont(FontProperties(family=family), fallback_to_default=False)
            return family
        except ValueError:
            pass
    return "DejaVu Sans"


def resolve_linear_layout(args):
    layout = getattr(args, "linear_layout", "legacy")
    publication = layout != "legacy" and getattr(args, "linear_style", "redwood") == "redwood"
    journal = getattr(args, "publication_journal", "nature")
    if args.width is None:
        widths = (88, 180) if journal == "nature-communications" else (89, 183)
        args.width = widths[layout == "two-column"] / 25.4 if publication else 13
    if args.max_reads is None:
        args.max_reads = 30 if publication and layout == "one-column" else 80
    if publication and args.width * 72 < 200:
        raise ValueError("Publication layouts require at least 200 pt width (2.78 inches)")
    return publication


def text_width(text, family, size=7, weight="normal"):
    return TextPath((0, 0), text, prop=FontProperties(family=family, size=size, weight=weight)).get_extents().width


def plan_feature_labels(features, length, width_pt, family, force_external=False):
    """Measure actual names, then allocate nonoverlapping external label lanes.

    The underlying feature intervals never move. Only text can be clamped to
    the plot boundary, with a leader connecting it to the correct feature.
    """
    result, lane_ends = [], []
    for feature in sorted(features, key=lambda f: (f["start"], f["stop"], f["name"])):
        start, stop = max(0, feature["start"]), min(length, feature["stop"])
        if stop <= start:
            continue
        head = min(max(length * .0045, 70), (stop - start) * .28, 240)
        available = (stop - start - head) / length * width_pt
        internal = not force_external and available >= text_width(feature["name"], family, weight="bold") + 1.5
        center = (start + stop + (head if feature["strand"] == "-" else -head)) / 2
        item = {"feature": feature, "internal": internal, "x": center, "lane": None}
        if not internal:
            size = text_width(feature["name"], family)
            if size > width_pt - 4:
                raise ValueError(f"Annotation name {feature['name']!r} exceeds the plot width; use a wider layout")
            target = (start + stop) / 2 / length * width_pt
            center_pt = min(width_pt - size / 2 - 1, max(size / 2 + 1, target))
            # A solitary short feature can use the empty space beside its
            # arrow, retaining the same physical label size without a new row.
            beside = stop / length * width_pt + 5 + size / 2
            if len(features) == 1 and beside + size / 2 < width_pt - 1:
                item.update(x=beside / width_pt * length, beside=True)
                result.append(item)
                continue
            left, right = center_pt - size / 2, center_pt + size / 2
            lane = len(lane_ends)
            for i, end in enumerate(lane_ends):
                shifted = max(center_pt, end + 4 + size / 2)
                if shifted + size / 2 <= width_pt - 1 and shifted - center_pt <= min(24, width_pt * .13):
                    lane, center_pt = i, shifted
                    left, right = center_pt - size / 2, center_pt + size / 2
                    break
            if lane == len(lane_ends):
                lane_ends.append(right)
            else:
                lane_ends[lane] = right
            item.update(x=center_pt / width_pt * length, lane=lane)
        result.append(item)
    return result, len(lane_ends)


def draw_publication_linear(args, reference, features, selected, evidence, rna):
    length, family = reference.length, publication_font()
    width_pt = args.width * 72
    gutter, right = 48, 7
    plot_width = width_pt - gutter - right
    fg, bg = ("#eef4fb", "#0d1117") if getattr(args, "dark", False) else ("#20252b", "white")
    colors = (BARK_COLOR, BARK_COLOR_ALT)
    requested = getattr(args, "linear_track", None)
    tracks = list(dict.fromkeys(requested if requested is not None else ["depth"]))
    if "none" in tracks:
        if len(tracks) > 1:
            raise ValueError("--linear-track none cannot be combined with other tracks")
        tracks = []
    if getattr(args, "hide_evidence", False):
        tracks = [track for track in tracks if track not in {"ends", "clips"}]
    class_mode = getattr(args, "read_color", "wood") == "class"
    if class_mode and not args.read_classes:
        raise ValueError("--read-color class requires --read-classes")
    palette = class_colors(evidence["depth"] if evidence else [])
    placed, row_count = pack_linear_reads(selected, length)
    nc = getattr(args, "publication_journal", "nature") == "nature-communications"
    thin, body, pitch = (1., 1.25, 2.) if nc else (.25, .70, 1.05)
    bands, plans = [], {}
    # Genome structure is read first; coverage stays immediately above reads.
    terminal = [f for f in features if f["type"] in {"repeat_region", "misc_feature"}]
    if terminal:
        bands.append(("termini", 14))
    trnas = [f for f in features if f["type"] == "tRNA"]
    if trnas:
        plans["trna"] = plan_feature_labels(trnas, length, plot_width, family, True)
        bands.append(("trna", 10 + 8 * plans["trna"][1]))
    genes = annotation_rows(features)
    for i, (_, row) in enumerate(genes):
        key = f"genes_{i}"
        plans[key] = plan_feature_labels(row, length, plot_width, family)
        bands.append((key, 10 + 8 * plans[key][1]))
    bands.append(("at", 17))
    if "gc" in getattr(args, "extra_tracks", []):
        bands.append(("gc", 17))
    if rna is not None:
        bands.append(("rna", 18))
    if evidence and "depth" in tracks:
        bands.append(("depth", 18))
    if row_count:
        bands.append(("reads", max(16, row_count * pitch + 4)))
    for track in ("ends", "clips"):
        if evidence and track in tracks:
            bands.append((track, 23))
    class_rows = int(np.ceil(len(palette) / (1 if args.linear_layout == "one-column" else 3))) if class_mode else 0
    footer = class_rows * 11 + (22 if getattr(args, "show_terminal_sequences", False) else 0)
    top_margin = 19
    height_pt = top_margin + sum(h for _, h in bands) + 4 + footer
    fig = plt.figure(figsize=(args.width, height_pt / 72))
    fig.patch.set_facecolor(bg)
    ax = fig.add_axes([gutter / width_pt, 0, plot_width / width_pt, 1])
    ax.set_xlim(.5, length + .5)
    ax.set_ylim(height_pt, 0)
    ax.set_axis_off()
    metadata = {"layout": args.linear_layout, "journal": args.publication_journal,
                "width_mm": args.width * 25.4, "height_mm": height_pt / 72 * 25.4,
                "font_family": family, "label_font_pt": 7, "secondary_font_pt": 6,
                "read_pitch_pt": pitch, "read_body_pt": body, "minimum_read_width_pt": thin,
                "read_rows": row_count, "tracks": tracks, "annotations": [], "composition": {}}

    def text(x, y, value, size=7, **kwargs):
        return ax.text(x, y, value, fontsize=size, family=family, color=kwargs.pop("color", fg),
                       va="center", **kwargs)

    def label(value, center, detail=None):
        x = -(7 / plot_width) * length
        text(x, center - (4 if detail else 0), value, ha="right", clip_on=False)
        if detail:
            text(x, center + 5, detail, size=6, ha="right", clip_on=False)

    # A single unit label and a sparse, shared, 1-based reference axis.
    label("Position", 9, "(kb)")
    target_ticks = 5 if args.linear_layout == "one-column" else 7
    raw = length / max(1, target_ticks - 1)
    power = 10 ** np.floor(np.log10(max(1, raw)))
    step = next(v * power for v in (1, 2, 2.5, 5, 10) if v * power >= raw)
    ticks = [1] + list(np.arange(step, length * .93, step)) + [length]
    end_width = text_width(f"{length / 1000:g}", family, 6.5)
    ticks = [pos for pos in ticks if pos in (1, length) or
             (length - pos) / length * plot_width > end_width + text_width(f"{pos / 1000:g}", family, 6.5) / 2 + 4]
    for pos in ticks:
        text(pos, 7, "1 bp" if pos == 1 else f"{pos / 1000:g}", size=6.5,
             ha="left" if pos == 1 else "right" if pos == length else "center")
        ax.vlines(pos, 12, 15, color=fg, linewidth=max(thin, .5))
    if getattr(args, "panel_label", None):
        fig.text(4 / width_pt, 1 - 8 / height_pt, args.panel_label,
                 fontsize=8 if not nc else 7, weight="bold", family=family, color=fg, va="center")

    top = top_margin
    for key, height in bands:
        if key == "termini":
            itr_caps = all((f["type"] == "repeat_region" and f["attributes"].get("rpt_type") == "inverted") or
                           (f["type"] == "misc_feature" and "cap" in f["name"].lower()) for f in terminal)
            label("ITR / cap" if itr_caps else "Features", top + height / 2)
            for f in sorted(terminal, key=lambda f: f["type"] != "repeat_region"):
                repeat = f["type"] == "repeat_region"
                add_arrow(ax, f["start"], f["stop"], top + 9, 4 if repeat else 7,
                          f["strand"], length, BARK_COLOR_ALT if repeat else REDWOOD_GRADIENT[0], 1)
            # Legend uses shapes and text, placed in the otherwise empty center.
            items = [("ITR" if itr_caps else "Repeat", BARK_COLOR_ALT, 4),
                     ("Cap" if itr_caps else "Feature", REDWOOD_GRADIENT[0], 7)]
            for offset, (name, color, h) in zip((-27, 9), items):
                x = length / 2 + offset / plot_width * length
                ax.add_patch(Rectangle((x, top + 4 - h / 2), 7 / plot_width * length, h,
                                       facecolor=color, edgecolor="none"))
                text(x + 10 / plot_width * length, top + 4, name, 6, ha="left")
        elif key in plans:
            plan, lanes = plans[key]
            center = top + lanes * 8 + 5
            label("tRNA" if key == "trna" else f"Genes ({genes[int(key.split('_')[1])][0]})", center)
            for item in plan:
                f = item["feature"]
                add_arrow(ax, f["start"], f["stop"], center, 8,
                          f["strand"], length, FEATURE_COLORS.get(f["type"], "#d08c35"), 1)
                if item["internal"]:
                    artist = text(item["x"], center, f["name"], weight="bold", color="white", ha="center")
                elif item.get("beside"):
                    artist = text(item["x"], center, f["name"], ha="center")
                else:
                    y = top + item["lane"] * 8 + 3
                    anchor = (max(0, f["start"]) + min(length, f["stop"])) / 2
                    ax.plot([item["x"], anchor, anchor], [y + 4, center - 6, center - 4],
                            color=fg, lw=max(thin, .5), solid_capstyle="butt")
                    artist = text(item["x"], y, f["name"], ha="center")
                artist.set_gid(f"annotation_label_{len(metadata['annotations'])}")
                metadata["annotations"].append({"name": f["name"], "type": f["type"],
                                                "start": f["start"] + 1, "end": f["stop"],
                                                "external": not item["internal"]})
        elif key in {"at", "gc"}:
            values = linear_base_fraction(reference.sequence, key.upper())
            lo, hi = np.percentile(values, [2, 98])
            scaled = (AT_RANGE[0] + np.clip((values - lo) / (hi - lo), 0, 1) * (AT_RANGE[1] - AT_RANGE[0])
                      if hi > lo else np.full_like(values, sum(AT_RANGE) / 2))
            ax.imshow(AT_COLORMAP(scaled)[None, :, :], extent=(.5, length + .5, top + 12, top + 5),
                      aspect="auto", interpolation="nearest")
            label(f"{key.upper()} (%)", top + 5)
            # Small key lives in the label gutter, without stealing genome width.
            left = -(40 / plot_width) * length
            end = -(7 / plot_width) * length
            gradient = AT_COLORMAP(np.linspace(*AT_RANGE, 128))[None, :, :]
            ax.imshow(gradient, extent=(left, end, top + 11, top + 8), aspect="auto", clip_on=False)
            text(left, top + 15, f"{lo * 100:.0f}", 6, ha="left", clip_on=False)
            text(end, top + 15, f"{hi * 100:.0f}", 6, ha="right", clip_on=False)
            metadata["composition"][key] = {"window_bp": 201, "limits_percent": [float(lo * 100), float(hi * 100)],
                                            "normalization_percentiles": [2, 98], "terminal_windows": "truncated"}
        elif key in {"rna", "depth"}:
            if key == "rna":
                forward, reverse = rna
                total, scale = forward + reverse, "log"
                rna_label = getattr(args, "rnaseq_label", "RNA depth")
                if text_width(rna_label, family) > gutter - 10:
                    # Keep long sample identifiers in the caption, where they
                    # remain editable and cannot squeeze the genome axis.
                    rna_label = "RNA depth"
                label(rna_label, top + 9, "log")
            else:
                total = sum(evidence["depth"].values(), np.zeros(length))
                scale = args.depth_scale
                label("Read depth", top + 9, "log" if scale == "log" else None)
            maximum = float(total.max()) if len(total) else 0
            transform = np.log1p if scale == "log" else lambda x: x
            heights = transform(total) / (float(transform(maximum)) or 1) * 10
            baseline = top + 17
            if key == "rna" and "rnaseq-strand" in getattr(args, "extra_tracks", []):
                ratio = np.divide(reverse, total, out=np.zeros_like(reverse, dtype=float), where=total > 0)
                rgba = np.array(to_rgba(colors[0], .85))[None, :] * (1 - ratio[:, None])
                rgba += np.array(to_rgba(colors[1], .85))[None, :] * ratio[:, None]
                polygons = [[(i + .5, baseline), (i + .5, baseline - h),
                             (i + 1.5, baseline - h), (i + 1.5, baseline)] for i, h in enumerate(heights)]
                ax.add_collection(PolyCollection(polygons, facecolors=rgba, edgecolors="none"))
            else:
                ax.fill_between(np.arange(length + 1) + .5, baseline, baseline - np.r_[heights, heights[-1]],
                                step="post", color=colors[0], alpha=.85, linewidth=0)
            text(length, top + 3, f"0-{maximum:,.0f}x", 6, ha="right")
            metadata[f"{key}_scale"] = {"transform": "log1p" if scale == "log" else "linear", "max_depth": maximum}
        elif key == "reads":
            label("Reads", top + height / 2, f"n = {len(placed)}")
            for _, group, row in placed:
                for segment in group:
                    gradient = [palette[segment.read_class]] if class_mode else REDWOOD_GRADIENT
                    add_gradient_read(ax, segment, top + 2 + (row + .5) * pitch, body,
                                      getattr(args, "min_indel", 10), gradient, min_width=thin)
        elif key in {"ends", "clips"}:
            profiles, clips = evidence["profiles"], key == "clips"
            first = profiles["short_left_clips"] + profiles["long_left_clips"] if clips else profiles["starts"]
            second = profiles["short_right_clips"] + profiles["long_right_clips"] if clips else profiles["ends"]
            x, above, widths = binned(first, args.bin_size)
            _, below, _ = binned(second, args.bin_size)
            maximum = max(int(above.max()), int(below.max()), 1)
            center, amplitude = top + 15, 6
            label("Soft clips" if clips else "Read ends", center)
            text(1, top + 5, "Above: left / below: right" if clips else "Above: starts / below: ends", 6, ha="left")
            text(length, top + 5, f"max {maximum:,}", 6, ha="right")
            for side, counts, sign, color in (("left", above, -1, colors[0]), ("right", below, 1, colors[1])):
                h = sign * np.log1p(counts) / np.log1p(maximum) * amplitude
                ax.bar(x, h, bottom=center, width=widths, align="edge", color=color, linewidth=0)
                if clips:
                    _, short, _ = binned(profiles[f"short_{side}_clips"], args.bin_size)
                    fraction = np.divide(short, counts, out=np.zeros_like(short, dtype=float), where=counts > 0)
                    ax.bar(x, h * fraction, bottom=center, width=widths, align="edge",
                           color=REDWOOD_GRADIENT[-1], linewidth=0)
            ax.hlines(center, .5, length + .5, color=colors[0], linewidth=thin, alpha=.45)
        top += height
    if class_mode:
        ncol = 1 if args.linear_layout == "one-column" else 3
        for i, (name, color) in enumerate(palette.items()):
            x = (i % ncol) * length / ncol
            y = top + 6 + (i // ncol) * 11
            ax.add_patch(Rectangle((x, y - 3), 6 / plot_width * length, 6, facecolor=color, edgecolor="none"))
            text(x + 9 / plot_width * length, y, name, 6, ha="left")
        top += class_rows * 11
    if getattr(args, "show_terminal_sequences", False):
        for y, sequence in ((top + 5, f"Left 5'-{reference.sequence[:30].upper()}..."),
                            (top + 15, f"Right ...{reference.sequence[-30:].upper()}-3'")):
            ax.text(1, y, sequence, family="monospace", fontsize=6, color=fg, va="center")
    # imshow can change limits even for the gutter color key.
    ax.set_xlim(.5, length + .5)
    ax.set_ylim(height_pt, 0)
    fig._redwood_publication_info = metadata
    if evidence:
        evidence["summary"]["production_style"] = {
            "tracks": tracks, "read_color": args.read_color, "read_rows": row_count,
            "at_window": 201, "composition_windows": "truncated at termini",
            "rnaseq_label": args.rnaseq_label if rna is not None else None,
            "publication": metadata,
        }
    return fig


def publication_caption(args, reference, evidence, info):
    title = getattr(args, "title", None) or reference.name
    lines = [f"{title}. Linear reference {reference.name}, {reference.length:,} bp.", "",
             "Gene arrow positions and directions follow the supplied annotation. Every visible gene and tRNA is labeled. "
             "ITRs and terminal caps, where shown, are supplied annotations, not inferred by Redwood."]
    if evidence:
        selection = evidence["summary"]["read_selection"]
        lines.extend(["", f"The read panel shows {len(selection['read_ids'])} unique reads selected by "
                      f"{selection['method']} ({selection['ranking']} ranking; {selection['terminal_window']} bp terminal windows), "
                      "then ordered by reference start and end. Supplementary segments stay with their read. "
                      f"Terminal groups: {selection['terminal_groups']}. The longest-read selection is intentional and is not a random population sample.",
                      "", f"Depth and endpoint/clipping profiles use all {evidence['summary']['reads']:,} unique input read names. "
                      "Depth counts aligned bases once per read, merging overlapping primary/supplementary intervals. "
                      "Endpoint and clipping profiles count alignment endpoints; starts are above and ends below the baseline. "
                      "Soft clips use reference-left above and reference-right below, independent of alignment strand. "
                      f"Light clipping segments are 1-{args.clip_threshold - 1} bp; dark segments are at least {args.clip_threshold} bp. "
                      f"Endpoint/clipping bins are {args.bin_size} bp, with height proportional to log(1 + count). "
                      "The displayed maximum is the largest binned count in either direction."])
    if args.read_color == "wood":
        lines += ["", "The brown read gradient is Redwood styling along increasing reference coordinates; it does not encode strand, quality or alleles. "
                  "Read-width changes mark CIGAR insertions and deletions, with a minimum printable width."]
    for key, values in info.get("composition", {}).items():
        lo, hi = values["limits_percent"]
        lines += ["", f"{key.upper()} composition uses centered 201 bp windows shortened at the termini. "
                  f"The color key spans the reference's 2nd-98th percentiles ({lo:.2f}-{hi:.2f}%); values outside that range are saturated."]
    if args.rnaseq_bam:
        lines += ["", f"{args.rnaseq_label}: primary RNA alignment base depth (not UMI counts), plotted as log(1 + depth). "
                  "Library preparation and input-population definitions should be described with the experiment."]
    if getattr(args, "subtitle", None):
        lines += ["", f"Sample information: {args.subtitle}."]
    lines += ["", f"Export: {info['width_mm']:.2f} x {info['height_mm']:.2f} mm; "
              f"{info['font_family']}; 7 pt essential labels, 6-6.5 pt secondary text. Place at 100% size."]
    return "\n".join(lines) + "\n"
