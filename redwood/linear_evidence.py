"""Coordinate-preserving inputs and descriptive evidence for linear genome plots.

All intervals here are zero-based, half-open. TSV positions and figure labels
are one-based. No statistic in this module is an automatic topology call.
"""

from __future__ import annotations

import csv
import gzip
import json
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pandas as pd
import pysam


@dataclass
class Reference:
    name: str
    sequence: str

    @property
    def length(self):
        return len(self.sequence)


def load_reference(path):
    """Require exactly one sequence, retaining case for masking diagnostics."""
    opener = gzip.open if str(path).endswith(".gz") else open
    names, chunks = [], []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                words = line[1:].split()
                if not words:
                    raise ValueError("FASTA sequence header must have a name")
                names.append(words[0])
            elif line:
                if not names:
                    raise ValueError("FASTA sequence must follow a named header")
                chunks.append(line)
    if len(names) != 1 or not chunks:
        raise ValueError("linear mode requires a FASTA with exactly one nonempty sequence")
    return Reference(names[0], "".join(chunks))


def load_classes(path):
    """Accept the issue's read/class/locus format and phase_reads_fast output."""
    classes = {}
    if path is None:
        return classes
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"read", "class"}.issubset(reader.fieldnames or []):
            raise ValueError("--read-classes needs a TSV header containing read and class")
        for row in reader:
            name, label = row["read"], row["class"]
            if not name or not label:
                raise ValueError("read/class values must not be empty")
            item = {"class": label, "locus": row.get("locus") or row.get("nuclear_locus") or ""}
            if name in classes and classes[name] != item:
                raise ValueError(f"conflicting read classes for {name!r}")
            classes[name] = item
    return classes


def load_annotations(path, reference):
    if path is None:
        return []
    features = []
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region"):
                fields = line.split()
                if len(fields) != 4 or fields[1:] != [reference.name, "1", str(reference.length)]:
                    raise ValueError(f"GFF sequence-region does not match FASTA: {line.strip()}")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GFF line {number}: expected nine tab-separated fields")
            name, _, kind, start, stop, _, strand, _, attrs = fields
            if name != reference.name:
                raise ValueError(f"GFF reference {name!r} does not match FASTA {reference.name!r}")
            start, stop = int(start) - 1, int(stop)
            if stop <= start:
                raise ValueError(f"GFF line {number}: end precedes start")
            attributes = dict(part.split("=", 1) for part in attrs.split(";") if "=" in part)
            label = next((attributes[k] for k in ("Name", "gene", "product", "ID") if attributes.get(k)), kind)
            features.append(dict(type=kind, start=start, stop=stop, strand=strand,
                                 name=unquote(label), attributes=attributes))
    # Many GFF writers repeat the same feature as gene + CDS/RNA, sometimes
    # with several redundant gene parents. Keep standalone genes (partial ITR
    # copies included), and prefer the specific type for an identical interval.
    specific = {(f["start"], f["stop"], f["strand"]) for f in features
                if f["type"] in {"CDS", "rRNA", "tRNA", "ncRNA"}}
    seen, result = set(), []
    for item in features:
        interval = (item["start"], item["stop"], item["strand"])
        if item["type"] in {"region", "source", "mRNA", "exon"}:
            continue
        if item["type"] == "gene" and interval in specific:
            continue
        key = (item["type"], *interval)
        if key not in seen:
            seen.add(key)
            result.append(item)
    return result


def merged_intervals(intervals):
    result = []
    for start, stop in sorted(intervals):
        if stop <= start:
            continue
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(stop, result[-1][1]))
        else:
            result.append((start, stop))
    return result


@dataclass
class Segment:
    name: str
    start: int
    stop: int
    reverse: bool
    supplementary: bool
    mapq: int
    cigar: list
    blocks: list
    query_start: int
    query_stop: int
    read_length: int
    left_clip: int
    right_clip: int
    hard_clipped: bool
    read_class: str
    locus: str


def read_segments(path, reference, classes=None):
    has_classes = classes is not None
    classes = classes or {}
    segments, counts = [], Counter(primary=0, supplementary=0, secondary=0, hard_clipped=0)
    with pysam.AlignmentFile(path, "rb") as bam:
        if reference.name not in bam.references:
            raise ValueError(f"BAM does not contain FASTA reference {reference.name!r}")
        length = bam.get_reference_length(reference.name)
        if length != reference.length:
            raise ValueError(f"BAM/FASTA length mismatch for {reference.name}: {length} != "
                             f"{reference.length}; linear mode requires an undoubled reference")
        if not bam.has_index():
            raise ValueError(f"BAM must be indexed: samtools index {path}")
        for read in bam.fetch(reference.name):
            if read.is_unmapped:
                continue
            if read.is_secondary:
                counts["secondary"] += 1
                continue
            counts["supplementary" if read.is_supplementary else "primary"] += 1
            cigar = read.cigartuples or []
            if not cigar or read.reference_end is None:
                raise ValueError(f"mapped read {read.query_name!r} has no reference-consuming CIGAR")
            if read.reference_start < 0 or read.reference_end > length:
                raise ValueError(f"alignment for {read.query_name!r} extends beyond the linear reference")
            hard = any(op == 5 for op, _ in cigar)
            counts["hard_clipped"] += int(hard)
            read_length = read.infer_read_length() or 0
            offset = cigar[0][1] if cigar[0][0] == 5 else 0
            qs, qe = read.query_alignment_start + offset, read.query_alignment_end + offset
            if read.is_reverse:
                qs, qe = read_length - qe, read_length - qs
            info = classes.get(read.query_name, {"class": "unclassified" if has_classes else "all", "locus": ""})
            segments.append(Segment(
                read.query_name, read.reference_start, read.reference_end, read.is_reverse,
                read.is_supplementary, read.mapping_quality, cigar, read.get_blocks(),
                qs, qe, read_length,
                cigar[0][1] if cigar[0][0] == 4 else 0,
                cigar[-1][1] if cigar[-1][0] == 4 else 0, hard,
                info["class"], info["locus"],
            ))
    return segments, dict(counts)


def representative_segment(group):
    """Use the longest primary alignment, or a supplementary one if necessary.

    Separate split segments must not masquerade as one end-to-end alignment.
    """
    return min(group, key=lambda s: (s.supplementary, -(s.stop - s.start), s.start, s.stop))


def terminal_group(segment, length, window=30):
    window = min(window, max(1, length // 2))
    left, right = segment.start < window, segment.stop > length - window
    return "both" if left and right else "left" if left else "right" if right else "internal"


def pack_linear_reads(segments, length):
    """Pack in coordinate order, keeping each read's split segments on one row."""
    groups = defaultdict(list)
    for segment in segments:
        groups[segment.name].append(segment)
    ordered = sorted(groups.items(), key=lambda item: (min(s.start for s in item[1]),
                                                       max(s.stop for s in item[1]), item[0]))
    rows, placed = [], []
    pad = length * .004
    for name, group in ordered:
        intervals = merged_intervals((max(0, s.start - pad), min(length, s.stop + pad)) for s in group)
        for row, occupied in enumerate(rows):
            if not any(s1 < e2 and s2 < e1 for s1, e1 in intervals for s2, e2 in occupied):
                occupied.extend(intervals)
                break
        else:
            row = len(rows)
            rows.append(list(intervals))
        placed.append((name, sorted(group, key=lambda s: (s.start, s.stop, s.supplementary)), row))
    return placed, len(rows)


def select_display_segments(segments, length, query=None, sort="ALNLEN", max_reads=80,
                            selection="terminal-balanced", terminal_window=30):
    """Select read IDs by length/terminus first, then order by reference position.

    Query clauses determine eligibility. All alignments of an eligible selected
    read are retained, including supplementary segments that fail the query.
    Display selection never filters the depth or evidence input.
    """
    if max_reads < 0:
        raise ValueError("--max-reads must be nonnegative (0 hides the read panel)")
    if terminal_window < 1:
        raise ValueError("terminal window must be positive")
    if selection not in {"terminal-balanced", "longest"}:
        raise ValueError(f"unknown linear read selection: {selection}")
    rows = []
    for segment in segments:
        cigar = segment.cigar
        rows.append({"POS": segment.start + 1, "MAPQ": segment.mapq,
                     "ALNLEN": sum(n for op, n in cigar if op not in {1, 4, 5}),
                     "MAPLEN": sum(n for op, n in cigar if op != 1),
                     "TRULEN": sum(n for _, n in cigar),
                     "READ": segment.name, "CLASS": segment.read_class, "reflength": length})
    if not rows or max_reads == 0:
        return []
    frame = pd.DataFrame(rows)
    # ['False'] was the original plotter's sentinel for disabling filters.
    if query and list(query) != ["False"]:
        for clause in query:
            frame = frame.query(clause, engine="python")
    groups = defaultdict(list)
    for segment in segments:
        groups[segment.name].append(segment)
    eligible = set(frame["READ"])
    representatives = {name: representative_segment(group) for name, group in groups.items() if name in eligible}
    # POS controls placement, never truncation. Selecting a coordinate prefix
    # would silently exclude reads anchored at the other end of the reference.
    rank = "ALNLEN" if sort == "POS" else sort
    excluded_ops = {"ALNLEN": {1, 4, 5}, "MAPLEN": {1}, "TRULEN": set()}[rank]
    ranked = sorted(representatives, key=lambda name: (
        -sum(n for op, n in representatives[name].cigar if op not in excluded_ops), name))
    if selection == "longest":
        names = set(ranked[:max_reads])
    else:
        buckets = {key: deque() for key in ("left", "right", "both", "internal")}
        for name in ranked:
            buckets[terminal_group(representatives[name], length, terminal_window)].append(name)
        names = set()
        # Exclusive groups count spanning reads once. Redistribute empty slots
        # among the remaining terminal groups before filling with internal reads.
        while len(names) < max_reads and any(buckets[key] for key in ("left", "right", "both")):
            for key in ("left", "right", "both"):
                if buckets[key] and len(names) < max_reads:
                    names.add(buckets[key].popleft())
        while len(names) < max_reads and buckets["internal"]:
            names.add(buckets["internal"].popleft())
    return sorted((s for s in segments if s.name in names),
                  key=lambda s: (s.start, s.stop, s.name, s.supplementary))


def classify_junction(left, right, length, window):
    gap = right.query_start - left.query_stop
    if left.hard_clipped or right.hard_clipped:
        return "unassessed_hard_clip"
    if abs(gap) > 500:
        return "unassessed_query_gap_or_overlap"
    # Use the query-facing breakpoints, not either end of either alignment.
    exit_pos = left.start if left.reverse else left.stop
    entry_pos = right.stop if right.reverse else right.start
    near_left = lambda pos: pos <= window
    near_right = lambda pos: pos >= length - window
    if left.reverse == right.reverse:
        wraps = (near_left(exit_pos) and near_right(entry_pos)) if left.reverse else (
            near_right(exit_pos) and near_left(entry_pos))
        return "end_to_start" if wraps else "internal_same_strand"
    terminal = ((near_left(exit_pos) and near_left(entry_pos)) or
                (near_right(exit_pos) and near_right(entry_pos)))
    return "terminal_inverted" if terminal else "internal_inverted"


def collect_evidence(segments, reference, counts, terminal_window=30, junction_window=300, clip_threshold=100):
    length = reference.length
    if min(terminal_window, junction_window, clip_threshold) < 1:
        raise ValueError("terminal/junction windows and clip threshold must be positive")
    terminal_window = min(terminal_window, max(1, length // 2))
    junction_window = min(junction_window, max(1, length // 2))
    groups = defaultdict(list)
    for segment in segments:
        groups[segment.name].append(segment)
    depth_diff = {}
    profile_names = ("starts", "ends", "short_left_clips", "long_left_clips",
                     "short_right_clips", "long_right_clips", "mapq0_starts", "mapq0_ends")
    profiles = {name: np.zeros(length, dtype=np.int64) for name in profile_names}
    class_reads, junctions, read_lengths = Counter(), [], []
    terminal_clips = {side: Counter(none=0, short=0, long=0, unassessed_hard_clip=0) for side in ("left", "right")}
    unaligned_reads = 0
    for name, group in groups.items():
        label = group[0].read_class
        class_reads[label] += 1
        diff = depth_diff.setdefault(label, np.zeros(length + 1, dtype=np.int64))
        # A molecule counts once at a reference base even where its primary and
        # supplementary alignments overlap. D/N gaps never contribute depth.
        for start, stop in merged_intervals(block for s in group for block in s.blocks):
            diff[start] += 1
            diff[stop] -= 1
        read_length = max(s.read_length for s in group)
        read_lengths.append(read_length)
        aligned_query = sum(e - s for s, e in merged_intervals((s.query_start, s.query_stop) for s in group))
        unaligned_reads += int(read_length - aligned_query >= 300)
        ordered = sorted(group, key=lambda s: (s.query_start, s.query_stop, s.start, s.reverse))
        for left, right in zip(ordered, ordered[1:]):
            junctions.append(dict(read=name, read_class=label, locus=left.locus,
                                  kind=classify_junction(left, right, length, junction_window),
                                  first_start=left.start + 1, first_end=left.stop,
                                  first_strand="-" if left.reverse else "+",
                                  second_start=right.start + 1, second_end=right.stop,
                                  second_strand="-" if right.reverse else "+",
                                  query_gap=right.query_start - left.query_stop))
        for s in group:
            profiles["starts"][s.start] += 1
            profiles["ends"][s.stop - 1] += 1
            if s.mapq == 0:
                profiles["mapq0_starts"][s.start] += 1
                profiles["mapq0_ends"][s.stop - 1] += 1
            for side, clip, pos, terminal in (
                ("left", s.left_clip, s.start, s.start < terminal_window),
                ("right", s.right_clip, s.stop - 1, s.stop > length - terminal_window),
            ):
                if clip:
                    category = "long" if clip >= clip_threshold else "short"
                    profiles[f"{category}_{side}_clips"][pos] += 1
                if terminal:
                    category = ("unassessed_hard_clip" if s.hard_clipped else
                                "none" if not clip else "long" if clip >= clip_threshold else "short")
                    terminal_clips[side][category] += 1
    depth = {label: np.cumsum(diff)[:-1] for label, diff in sorted(depth_diff.items())}
    total = sum(depth.values(), np.zeros(length, dtype=np.int64))
    middle = total[length // 4:max(length // 4 + 1, 3 * length // 4)]
    shoulder = max(1, int(length * .02))
    mid_mean = float(np.mean(middle))
    n = len(groups)
    alignments = len(segments)
    warnings = []
    if not n:
        warnings.append("No mapped primary/supplementary reads on the reference.")
    if counts["hard_clipped"]:
        warnings.append("Hard-clipped alignments present: clip/junction evidence is incomplete; remap with minimap2 -Y.")
    if not counts["supplementary"]:
        warnings.append("No supplementary alignments present: missing split joins cannot establish their absence in the reads.")
    if counts["secondary"]:
        warnings.append("Secondary alignments excluded from all tracks and evidence.")
    masked = sum(base.islower() for base in reference.sequence) / length
    if masked:
        warnings.append("Reference contains lowercase (soft-masked) bases.")
    if n and unaligned_reads / n > .05:
        warnings.append("More than 5% of reads have >=300 bp outside mitochondrial alignments; NUMTs/chimeras require nuclear co-mapping to assess.")
    if any(label != "all" for label in class_reads):
        warnings.append("Read classes are supplied labels; identical recent NUMTs cannot be separated from mtDNA by SNPs alone.")
    summary = {
        "reference": reference.name, "length": length, "requested_topology": "linear",
        "topology_call": "not_assessed", "warnings": warnings,
        "alignment_counts": counts, "reads": n, "lowercase_fraction": masked,
        "mapq0_alignment_fraction": sum(s.mapq == 0 for s in segments) / alignments if alignments else None,
        "reads_with_300bp_unaligned": unaligned_reads,
        "unaligned_read_fraction": unaligned_reads / n if n else None,
        "terminal_window": terminal_window, "junction_window": junction_window,
        "clip_threshold": clip_threshold,
        "terminal_starts": int(profiles["starts"][:terminal_window].sum()),
        "terminal_ends": int(profiles["ends"][-terminal_window:].sum()),
        "uniform_endpoint_expectation": alignments * terminal_window / length,
        "terminal_soft_clips": terminal_clips,
        "junction_counts": dict(Counter(row["kind"] for row in junctions)),
        "read_length_median": float(np.median(read_lengths)) if n else None,
        "read_length_p99": float(np.percentile(read_lengths, 99)) if n else None,
        "reads_longer_than_reference": sum(value > length for value in read_lengths),
        "middle_mean_depth": mid_mean,
        "left_shoulder_to_middle": float(np.mean(total[:shoulder])) / mid_mean if mid_mean else None,
        "right_shoulder_to_middle": float(np.mean(total[-shoulder:])) / mid_mean if mid_mean else None,
        "classes": {label: {"reads": class_reads[label], "mean_depth": float(d.mean())}
                    for label, d in depth.items()},
        "depth_definition": "unique reads per aligned reference base (M, =, X); overlapping segments merged; D/N excluded",
        "endpoint_definition": "alignment starts/ends in reference coordinates, including supplementary segments",
    }
    return dict(summary=summary, depth=depth, profiles=profiles, junctions=junctions)


def write_evidence(base, evidence):
    Path(f"{base}.evidence.json").write_text(json.dumps(evidence["summary"], indent=2) + "\n")
    columns = {**{f"depth:{label}": values for label, values in evidence["depth"].items()}, **evidence["profiles"]}
    with open(f"{base}.profiles.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["position", *columns])
        writer.writerows((i + 1, *(int(values[i]) for values in columns.values()))
                         for i in range(evidence["summary"]["length"]))
    fields = ["read", "read_class", "locus", "kind", "first_start", "first_end", "first_strand",
              "second_start", "second_end", "second_strand", "query_gap"]
    with open(f"{base}.junctions.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(evidence["junctions"])
