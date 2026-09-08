"""Optional selected-site allele companion for a linear genome figure.

Sites are supplied by the user, not discovered here. Counts describe the input
BAM population and are not variant validation, phasing, or heteroplasmy calls.
"""

from __future__ import annotations

from bisect import bisect_left
from collections import Counter, defaultdict
import csv
import math

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Rectangle

from .linear_evidence import pack_linear_reads, select_display_segments


BASE_COLORS = {"A": "#2e7d32", "C": "#1565c0", "G": "#a77700", "T": "#c62828"}
STATUSES = ("A", "C", "G", "T", "deletion", "conflict", "low_quality",
            "ambiguous", "missing_sequence", "skip", "uncovered")


def load_variant_sites(path, reference):
    """Read one-based position[/label] TSV, or a named pos/position column.

The supplied Hydra table is accepted directly. Optional reference/ref/
v5_allele(hap1) columns are checked against the actual FASTA; reported source
frequencies and inferred haplotype labels are deliberately not imported.
"""
    with open(path, newline="") as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t")
                if row and row[0].strip() and not row[0].lstrip().startswith("#")]
    if not rows:
        raise ValueError("variant sites file is empty")
    header = [value.strip().lower() for value in rows[0]]
    named = any(value in {"pos", "position"} for value in header)
    if named:
        position_index = next(i for i, value in enumerate(header) if value in {"pos", "position"})
        label_index = next((header.index(key) for key in ("label", "gene") if key in header), None)
        reference_index = next((header.index(key) for key in ("reference", "ref", "v5_allele(hap1)")
                                if key in header), None)
        rows = rows[1:]
    else:
        position_index, label_index, reference_index = 0, 1, None
    sites, seen = [], set()
    for row in rows:
        try:
            position = int(row[position_index])
        except (ValueError, IndexError) as exc:
            raise ValueError("variant sites need integer one-based pos/position values") from exc
        if not 1 <= position <= reference.length:
            raise ValueError(f"variant position {position} is outside reference 1–{reference.length}")
        if position in seen:
            raise ValueError(f"duplicate variant position {position}")
        seen.add(position)
        base = reference.sequence[position - 1].upper()
        if reference_index is not None:
            supplied = row[reference_index].strip().upper() if reference_index < len(row) else ""
            if supplied != base:
                raise ValueError(f"variant reference allele at {position}: {supplied!r} does not match FASTA {base}")
        label = row[label_index].strip() if label_index is not None and label_index < len(row) else ""
        sites.append({"position": position, "reference": base, "label": label})
    if not sites:
        raise ValueError("variant sites file contains no sites")
    return sorted(sites, key=lambda site: site["position"])


def _alignment_calls(read, positions, minimum_quality):
    """Inspect aligned query/reference pairs; stored BAM bases face the reference."""
    first = bisect_left(positions, read.reference_start)
    last = bisect_left(positions, read.reference_end)
    wanted = set(positions[first:last])
    if not wanted:
        return {}
    # Older supported pysam releases lack get_aligned_pairs(with_cigar=True).
    # Track N operations separately so a skipped region is never a deletion.
    skips, cursor = set(), read.reference_start
    for operation, span in read.cigartuples or []:
        if operation == 3:
            skips.update(position for position in wanted if cursor <= position < cursor + span)
        if operation in {0, 2, 3, 7, 8}:
            cursor += span
    sequence, qualities = read.query_sequence, read.query_qualities
    calls = {}
    for query_position, reference_position in read.get_aligned_pairs(matches_only=False):
        if reference_position not in wanted:
            continue
        if reference_position in skips:
            status, quality = "skip", None
        elif query_position is None:
            status, quality = "deletion", None
        elif sequence is None:
            status, quality = "missing_sequence", None
        else:
            base = sequence[query_position].upper()
            quality = int(qualities[query_position]) if qualities is not None else None
            status = ("low_quality" if quality is None or quality < minimum_quality else
                      base if base in BASE_COLORS else "ambiguous")
        calls[reference_position] = {"status": status, "base_quality": quality,
                                     "supplementary": bool(read.is_supplementary)}
    return calls


def _resolve_calls(observations):
    if not observations:
        return "uncovered"
    confident = [call for call in observations if call["status"] in {*BASE_COLORS, "deletion"}]
    if len({call["status"] for call in confident}) > 1:
        return "conflict"
    # Agreeing overlap counts once. A supplementary call can fill a missing or
    # low-quality primary call; it never silently overwrites a disagreement.
    candidates = confident or observations
    return min(candidates, key=lambda call: (call["supplementary"],
                                            -(call["base_quality"] if call["base_quality"] is not None else -1)))['status']


def collect_variant_calls(path, reference, sites, minimum_quality=20):
    """Return per-read/site calls and all-input-population counts without pileup caps."""
    import pysam

    if not isinstance(minimum_quality, int) or not 0 <= minimum_quality <= 93:
        raise ValueError("variant minimum base quality must be an integer from 0 to 93")
    positions = [site["position"] - 1 for site in sites]
    observations = defaultdict(lambda: defaultdict(list))
    names, mapq_zero_names = set(), set()
    alignment_counts = Counter(primary=0, supplementary=0, secondary_excluded=0)
    with pysam.AlignmentFile(path, "rb") as bam:
        if reference.name not in bam.references:
            raise ValueError(f"variant BAM does not contain FASTA reference {reference.name!r}")
        if bam.get_reference_length(reference.name) != reference.length:
            raise ValueError("variant BAM/FASTA length mismatch; use the undoubled reference")
        if not bam.has_index():
            raise ValueError(f"variant BAM must be indexed: samtools index {path}")
        for read in bam.fetch(reference.name):
            if read.is_unmapped:
                continue
            if read.is_secondary:
                alignment_counts["secondary_excluded"] += 1
                continue
            if read.reference_end is None or not read.cigartuples:
                raise ValueError(f"mapped read {read.query_name!r} lacks a reference-consuming CIGAR")
            names.add(read.query_name)
            if read.mapping_quality == 0:
                mapq_zero_names.add(read.query_name)
            alignment_counts["supplementary" if read.is_supplementary else "primary"] += 1
            for position, call in _alignment_calls(read, positions, minimum_quality).items():
                observations[read.query_name][position].append(call)
    calls = {name: [_resolve_calls(observations[name][position]) for position in positions]
             for name in sorted(names)}
    summaries = []
    for column, site in enumerate(sites):
        counts = Counter(row[column] for row in calls.values())
        called = sum(counts[base] for base in BASE_COLORS)
        summaries.append({**site, "counts": {status: counts[status] for status in STATUSES},
                          "called_reads": called,
                          "allele_fractions": {base: counts[base] / called if called else None
                                               for base in BASE_COLORS}})
    return calls, {"sites": summaries, "population_reads": len(names),
                   "reads_with_mapq_zero_alignment": len(mapq_zero_names),
                   "alignment_counts": dict(alignment_counts), "minimum_base_quality": minimum_quality,
                   "fraction_denominator": "unique reads with an unambiguous quality-filtered A/C/G/T call at the site",
                   "overlap_policy": "count each read once; prefer primary, fill unavailable primary from supplementary; conflicting confident calls are excluded",
                   "population_policy": "all mapped nonsecondary reads on the supplied reference; no additional class or MAPQ filtering",
                   "interpretation": "observed alleles at supplied sites; no variant validation, inferred haplotypes, or heteroplasmy call"}


def draw_variant_details(args, reference, selected):
    """Draw an optional native-size allele matrix and attach numerical metadata."""
    if not getattr(args, "variant_sites", None):
        raise ValueError("variant companion needs a --variant-sites TSV")
    if not getattr(args, "main_bam", None):
        raise ValueError("variant companion needs a long-read BAM")
    sites = load_variant_sites(args.variant_sites, reference)
    width = float(args.width)
    # A seven-site matrix fits at a one-column publication width without
    # squeezing essential 7 pt letters or coordinate labels.
    if not math.isfinite(width) or width * 72 - 38 < len(sites) * 24:
        needed = (len(sites) * 24 + 38) / 72
        raise ValueError(f"{len(sites)} variant columns need at least {needed:.2f} inches; increase --width or supply fewer sites")
    minimum_quality = getattr(args, "variant_min_base_quality", 20)
    calls, metadata = collect_variant_calls(args.main_bam, reference, sites, minimum_quality)
    limited = select_display_segments(selected, reference.length, max_reads=30,
                                      selection=getattr(args, "linear_read_selection", "terminal-balanced"),
                                      terminal_window=getattr(args, "terminal_window", 30))
    display_names = [name for name, _, _ in pack_linear_reads(limited, reference.length)[0]]
    selected_names = [name for name, _, _ in pack_linear_reads(selected, reference.length)[0]]
    absent = set(selected_names) - calls.keys()
    if absent:
        raise ValueError("selected variant-panel read IDs are absent from the input BAM")
    metadata.update({"reference": {"name": reference.name, "length": reference.length},
                     "source_bam": str(args.main_bam), "source_sites": str(args.variant_sites),
                     "display_policy": "up to 30 longest terminal-balanced reads from the selected pool, then coordinate order"
                     if getattr(args, "linear_read_selection", "terminal-balanced") == "terminal-balanced" else
                     "up to 30 longest reads from the selected pool, then coordinate order",
                     "selected_reads": [{"read": name, "calls": calls[name]} for name in selected_names],
                     "displayed_reads": [{"index": index + 1, "read": name, "calls": calls[name]}
                                         for index, name in enumerate(display_names)],
                     "selected_count": len(selected_names), "displayed_count": len(display_names),
                     "display_symbols": {"A": "A", "C": "C", "G": "G", "T": "T", "deletion": "−",
                                         "conflict": "!", "other": "·"}})
    dark = getattr(args, "dark", False)
    background, foreground = ("#0d1117", "#eef4fb") if dark else ("white", "#111827")
    muted = "#b6bfcc" if dark else "#475467"
    # Coordinates are in physical points. Each displayed read occupies 9 pt,
    # independent of output width, so a native-size export preserves its type.
    row_height, matrix_top = 9, 147
    height_points = matrix_top + max(1, len(display_names)) * row_height + 41
    fig = plt.figure(figsize=(width, height_points / 72), facecolor=background)
    ax = fig.add_axes([28 / (width * 72), 0, 1 - 38 / (width * 72), 1], facecolor=background)
    ax.set_xlim(0, len(sites))
    ax.set_ylim(height_points, 0)
    ax.set_axis_off()
    text = dict(color=foreground, fontsize=7, va="center")
    ax.text(0, 11, "Alleles at supplied sites", fontweight="bold", **text)
    ax.text(0, 23, f"All input reads: {metadata['population_reads']:,} · BQ ≥{minimum_quality}",
            color=muted, fontsize=6, va="center")
    for column, site in enumerate(metadata["sites"]):
        x = column + .5
        ax.text(x, 39, f"{site['position']:,}", ha="center", **text)
        label = site["label"]
        max_characters = max(4, int((width * 72 - 38) / len(sites) / 3.5))
        if len(label) > max_characters:
            label = label[:max_characters - 1] + "…"
        ax.text(x, 50, label, ha="center", color=muted, fontsize=6, va="center")
        ax.text(x, 62, f"Ref {site['reference']}", ha="center", **text)
        for index, base in enumerate(BASE_COLORS):
            y = 77 + index * 13
            fraction = site["allele_fractions"][base]
            value = f"{site['counts'][base]:,}\n{100 * fraction:.1f}%" if fraction is not None else "0\n—"
            ax.text(x, y, value, ha="center", fontsize=6, color=foreground, va="center", linespacing=1.0)
    for index, base in enumerate(BASE_COLORS):
        ax.text(-.10, 77 + index * 13, base, ha="right", fontweight="bold", **text)
    ax.text(0, 133, f"Selected reads: {len(display_names)} / {len(selected_names)} shown", **text)
    for index, name in enumerate(display_names):
        y = matrix_top + index * row_height
        ax.text(-.10, y + row_height / 2, str(index + 1), ha="right", color=muted, fontsize=6, va="center")
        for column, status in enumerate(calls[name]):
            if status in BASE_COLORS:
                rgb = to_rgb(BASE_COLORS[status])
                fill = tuple(.22 * value + .78 for value in rgb)
                symbol = status
            else:
                fill = "#dfe3e8" if status == "conflict" else "#f4f5f7"
                symbol = "−" if status == "deletion" else "!" if status == "conflict" else "·"
            ax.add_patch(Rectangle((column + .025, y + .35), .95, row_height - .7,
                                   facecolor=fill, edgecolor="none"))
            ax.text(column + .5, y + row_height / 2, symbol, ha="center", color="#111827", fontsize=7, va="center")
    if not display_names:
        ax.text(0, matrix_top + row_height / 2, "No selected reads", **text)
    bottom = matrix_top + max(1, len(display_names)) * row_height
    ax.text(0, bottom + 11, "Fractions: callable A/C/G/T reads at each site.", color=muted, fontsize=6, va="center")
    ax.text(0, bottom + 22, "· no call   − deletion   ! conflicting alignments", color=muted, fontsize=6, va="center")
    ax.text(0, bottom + 33, "Columns are equally spaced; read IDs are in the metadata.", color=muted, fontsize=6, va="center")
    fig._redwood_variant_metadata = metadata
    fig._redwood_caption = (
        f"Observed alleles at {len(sites)} supplied sites on {reference.name}. "
        f"Counts use all {metadata['population_reads']:,} unique nonsecondary input read names; "
        f"A/C/G/T calls require base quality at least {minimum_quality}. "
        "Overlapping primary and supplementary alignments count once per read. "
        "Conflicting confident calls are excluded from allele fractions, whose denominator is "
        "the number of callable A/C/G/T reads at each site. "
        f"The matrix shows {len(display_names)} of {len(selected_names)} selected reads, "
        "retaining the main plot's length-ranking/terminal-selection rule and coordinate ordering. "
        "Column spacing is categorical. Source-site frequencies and haplotype labels are not imported. "
        "Read IDs, detailed missing-call categories and numeric counts are recorded in the companion JSON and TSV files.\n"
    )
    return fig
