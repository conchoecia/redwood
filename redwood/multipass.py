"""Classification and layout of circular-genome reads.

A read that traverses the circular reference more than once in a single
continuous alignment is a candidate rolling-circle-amplification (RCA)
molecule. These multi-pass reads are drawn as inward spirals; every other
read is a regular single-pass read.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pysam


# Distinct shades of red for successive multi-pass reads (cycled).
MULTIPASS_COLORS = [
    "#d1322b", "#8c1d40", "#e8651e", "#a82339",
    "#f0883e", "#6f1d1b", "#cf4b6b", "#b5341d",
]

# CIGAR operations (pysam codes) that introduce a gap in the alignment.
_GAP_OPS = (1, 2, 3)  # insertion, deletion, reference skip


@dataclass
class MultiPassRead:
    """A read that circles the reference >1x in one continuous alignment."""
    start: int       # 0-based start on the single-copy circle
    passes: float    # number of circle traversals (>1)
    span: int        # reference span in bp
    # pysam CIGAR tuples in (op_code, length) order — for indel rendering.
    cigar: list[tuple[int, int]] = field(default_factory=list)
    marks: list[tuple[int, str]] = field(default_factory=list)   # (offset along the reference from start, kind)
    name: str = ""


@dataclass
class RegularRead:
    """A single-pass read."""
    start: int       # 0-based start on the single-copy circle
    span: int        # reference span, capped to one circle for drawing
    # pysam CIGAR tuples in (op_code, length) order — kept so the renderer can
    # draw per-read insertions and deletions.
    cigar: list[tuple[int, int]] = field(default_factory=list)
    marks: list[tuple[int, str]] = field(default_factory=list)   # (offset along the reference from start, kind)
    name: str = ""


def read_marks(read, reference_seq: str, true_length: int, min_indel: int = 1) -> list[tuple[int, str]]:
    """IGV-style marks of one alignment against the reference: ``(offset, kind)`` with ``offset`` the
    reference distance from the alignment start and ``kind`` the read base at a mismatch (A/C/G/T),
    ``I`` for an insertion of >= ``min_indel`` bp (placed at the preceding reference base) or ``D``
    for a deletion of >= ``min_indel`` bp (placed at its first deleted base). Positions are folded
    onto the single-copy circle so doubled references work."""
    marks: list[tuple[int, str]] = []
    seq = read.query_sequence
    if not seq:
        return marks
    start = read.reference_start
    ins_run = 0
    last_ref = None
    del_run = 0
    del_start = None
    for qpos, rpos in read.get_aligned_pairs():
        if rpos is None:            # insertion (or soft clip: qpos with no reference)
            if last_ref is not None:
                ins_run += 1
            continue
        if ins_run:
            if ins_run >= min_indel:
                marks.append((last_ref - start, "I"))
            ins_run = 0
        if qpos is None:            # deletion / reference skip
            if del_run == 0:
                del_start = rpos
            del_run += 1
            last_ref = rpos
            continue
        if del_run:
            if del_run >= min_indel:
                marks.append((del_start - start, "D"))
            del_run = 0
        base = seq[qpos].upper()
        if base in "ACGT" and base != reference_seq[rpos % true_length].upper():
            marks.append((rpos - start, base))
        last_ref = rpos
    if del_run >= min_indel and del_start is not None:
        marks.append((del_start - start, "D"))
    return marks


def classify_circular_reads(
    bam_path: Path | None,
    true_length: int,
    *,
    max_internal_gap: int = 50,
    min_pass_fraction: float = 1.0,
    reference_seq: str | None = None,
    min_indel: int = 1,
) -> tuple[list[MultiPassRead], list[RegularRead]]:
    """Split primary alignments into multi-pass (RCA) and regular reads.

    A read is multi-pass when (1) its alignment is continuous — no insertion,
    deletion, or reference skip longer than ``max_internal_gap`` bp, since a
    gap that large is treated as un-trimmed adapter or chimeric junk — and
    (2) its reference span covers at least ``min_pass_fraction`` circles.
    """
    multipass: list[MultiPassRead] = []
    regular: list[RegularRead] = []
    if bam_path is None or not Path(bam_path).exists():
        return multipass, regular
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            span = read.reference_length
            if not span:
                continue
            start = read.reference_start % true_length
            cigar = list(read.cigartuples or [])
            biggest_gap = 0
            for op, oplen in cigar:
                if op in _GAP_OPS and oplen > biggest_gap:
                    biggest_gap = oplen
            continuous = biggest_gap <= max_internal_gap
            marks = read_marks(read, reference_seq, true_length, min_indel) if reference_seq else []
            if continuous and span >= min_pass_fraction * true_length:
                multipass.append(
                    MultiPassRead(start, span / true_length, span, cigar, marks, read.query_name))
            else:
                regular.append(RegularRead(start, min(span, true_length), cigar, marks, read.query_name))
    # draw the longest spirals first (outermost)
    multipass.sort(key=lambda r: r.passes, reverse=True)
    return multipass, regular


def _arc_intervals(start: float, span: float, length: int) -> list[tuple[float, float]]:
    """Linear [lo, hi) sub-intervals in [0, length) covered by a circular arc."""
    if span >= length:
        return [(0.0, float(length))]
    s = start % length
    e = s + span
    if e <= length:
        return [(s, e)]
    return [(s, float(length)), (0.0, e - length)]


def _intervals_overlap(a, b) -> bool:
    return any(lo1 < hi2 and lo2 < hi1 for lo1, hi1 in a for lo2, hi2 in b)


def pack_circular_reads(reads, length: int, pad: float = 0.0):
    """Greedily pack regular reads onto rungs, longest first.

    Reads are placed in descending span order; each read goes on the lowest
    (outermost) rung where it does not overlap an already-placed read. So the
    longest unplaced read opens each new rung — outer rungs carry the longest
    reads and shorter reads backfill the gaps — rather than producing a
    start-coordinate cascade. Reads of equal span are ordered by start position,
    so a set of full-length reads (each on its own rung) forms one clean cascade
    around the circle. Returns ``(placed, n_rungs)`` where ``placed`` is a list
    of ``(RegularRead, rung_index)``.
    """
    rungs: list[list] = []          # rungs[i] = list of placed arc-interval sets
    placed: list[tuple] = []
    for rd in sorted(reads, key=lambda r: (-r.span, r.start)):
        ivals = _arc_intervals(rd.start - pad, rd.span + 2 * pad, length)
        for ri, occupied in enumerate(rungs):
            if not any(_intervals_overlap(ivals, o) for o in occupied):
                occupied.append(ivals)
                placed.append((rd, ri))
                break
        else:
            rungs.append([ivals])
            placed.append((rd, len(rungs) - 1))
    return placed, len(rungs)
