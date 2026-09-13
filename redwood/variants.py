"""Per-column mismatch / insertion / deletion table of mapped reads against the reference.

Every reference column gets the counts of A/C/G/T, deletions (reads with a gap at the column) and
insertions (reads carrying extra bases after the column). Positions are folded modulo the reference
length, so a BAM mapped to the doubled / tripled reference that ``redwood run`` builds for circular
genomes reports one row per position of the single-copy molecule. This is a description of the read
population, not a variant caller: a column is flagged when the read majority disagrees with the
reference (``mismatch`` / ``deletion``), when a fraction of reads carry an insertion (``insertion``),
or when a minor allele reaches ``min_minor_frac`` (``minor``, e.g. heteroplasmy, mixed haplotypes,
NUMT-derived reads, or a systematic sequencing error).
"""

from __future__ import annotations

import collections
import json
from pathlib import Path

import pysam

BASES = ("A", "C", "G", "T")
COLUMNS = [
    "pos", "ref", "depth", "A", "C", "G", "T", "del", "ins", "ins_seq",
    "major", "major_frac", "minor_frac", "ins_frac", "events",
]


def read_first_sequence(path: Path) -> tuple[str, str]:
    name, chunks = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                break
            name = line[1:].split()[0]
            continue
        if name is not None:
            chunks.append(line.strip())
    if name is None:
        raise ValueError(f"no FASTA record in {path}")
    return name, "".join(chunks).upper()


def _keep(read) -> bool:
    return not (read.is_unmapped or read.is_secondary or read.is_qcfail or read.is_duplicate)


def column_variants(
    bam_path: Path,
    reference_seq: str,
    *,
    contig: str | None = None,
    min_base_quality: int = 20,
    min_depth: int = 5,
    min_minor_frac: float = 0.05,
    fold: bool = True,
    max_depth: int = 1_000_000,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Return (rows, summary). ``rows`` has one dict per reference column (all columns, flagged or not).

    Base counts come from ``pysam.count_coverage`` (C speed, base quality >= ``min_base_quality``);
    deletions and insertions from one pass over the CIGAR strings, so the table is fast even at
    several thousand-fold depth."""
    length = len(reference_seq)
    if length == 0:
        raise ValueError("empty reference sequence")
    base_counts = {b: [0] * length for b in BASES}
    del_counts = [0] * length
    insertions = [collections.Counter() for _ in range(length)]
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        contigs = [contig] if contig else list(bam.references)
        for name in contigs:
            ref_len = bam.get_reference_length(name)
            if not fold and ref_len != length:
                continue
            if fold and ref_len % length != 0:
                raise ValueError(
                    f"BAM contig {name} is {ref_len} bp, not a multiple of the reference length {length}; "
                    "map to the single-copy or the doubled reference"
                )
            cov = bam.count_coverage(name, 0, ref_len, quality_threshold=min_base_quality, read_callback=_keep)
            for b, arr in zip(BASES, cov):
                target = base_counts[b]
                for i, v in enumerate(arr):
                    if v:
                        target[i % length] += int(v)
            for read in bam.fetch(name, 0, ref_len):
                if not _keep(read) or read.cigartuples is None:
                    continue
                seq = read.query_sequence
                rpos, qpos, last_ref = read.reference_start, 0, None
                for op, n in read.cigartuples:
                    if op in (0, 7, 8):          # M = X
                        rpos += n; qpos += n; last_ref = rpos - 1
                    elif op == 2:                # D
                        for k in range(n):
                            del_counts[(rpos + k) % length] += 1
                        rpos += n; last_ref = rpos - 1
                    elif op == 3:                # N
                        rpos += n; last_ref = rpos - 1
                    elif op == 1:                # I
                        if last_ref is not None and seq:
                            insertions[last_ref % length][seq[qpos:qpos + n].upper()] += 1
                        qpos += n
                    elif op == 4:                # S
                        qpos += n
    rows = []
    for pos in range(length):
        c = {b: base_counts[b][pos] for b in BASES}
        c["del"] = del_counts[pos]
        depth = sum(c.values())
        ref = reference_seq[pos]
        ins_total = sum(insertions[pos].values())
        ins_seq = insertions[pos].most_common(1)[0][0] if ins_total else ""
        if depth == 0:
            rows.append({"pos": pos + 1, "ref": ref, "depth": 0, "A": 0, "C": 0, "G": 0, "T": 0, "del": 0,
                         "ins": 0, "ins_seq": "", "major": "", "major_frac": 0.0, "minor_frac": 0.0,
                         "ins_frac": 0.0, "events": "nodepth"})
            continue
        major, major_n = max(((b, c[b]) for b in BASES + ("del",)), key=lambda item: (item[1], item[0] == ref))
        major_frac = major_n / depth
        ins_frac = ins_total / depth
        events = []
        if depth >= min_depth:
            if major == "del":
                events.append("deletion")
            elif major != ref:
                events.append("mismatch")
            if ins_frac >= min_minor_frac:
                events.append("insertion")
            if (1 - major_frac) >= min_minor_frac:
                events.append("minor")
        else:
            events.append("lowdepth")
        rows.append({
            "pos": pos + 1, "ref": ref, "depth": depth, "A": c["A"], "C": c["C"], "G": c["G"], "T": c["T"],
            "del": c["del"], "ins": ins_total, "ins_seq": ins_seq, "major": major,
            "major_frac": round(major_frac, 4), "minor_frac": round(1 - major_frac, 4),
            "ins_frac": round(ins_frac, 4), "events": ",".join(events),
        })
    return rows, summarize_variants(rows, length, min_minor_frac, min_depth)


def summarize_variants(rows: list[dict[str, object]], length: int, min_minor_frac: float, min_depth: int) -> dict[str, object]:
    covered = [r for r in rows if int(r["depth"]) > 0]
    def has(tag): return sum(1 for r in rows if tag in str(r["events"]).split(","))
    return {
        "columns": length,
        "columns_with_depth": len(covered),
        "columns_below_min_depth": has("lowdepth"),
        "mean_depth": round(sum(int(r["depth"]) for r in rows) / length, 2) if length else 0,
        "mismatch_columns": has("mismatch"),
        "deletion_columns": has("deletion"),
        "insertion_columns": has("insertion"),
        "minor_allele_columns": has("minor"),
        "min_minor_frac": min_minor_frac,
        "min_depth": min_depth,
    }


def flagged_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [r for r in rows if r["events"] and r["events"] not in ("nodepth", "lowdepth")]


def write_variant_table(rows: list[dict[str, object]], path: Path, all_columns: bool = False) -> int:
    selected = rows if all_columns else flagged_rows(rows)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in selected:
            handle.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")
    return len(selected)


def write_variants(
    bam_path: Path,
    mito_fasta: Path,
    output: Path,
    *,
    summary_path: Path | None = None,
    all_columns: bool = False,
    **kwargs,
) -> dict[str, object]:
    _, sequence = read_first_sequence(Path(mito_fasta))
    rows, summary = column_variants(Path(bam_path), sequence, **kwargs)
    summary["rows_written"] = write_variant_table(rows, Path(output), all_columns=all_columns)
    summary["table"] = str(output)
    if summary_path:
        Path(summary_path).write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def run_variants(args) -> dict[str, object]:
    summary = write_variants(
        Path(args.bam), Path(args.mito_fasta), Path(args.output), summary_path=args.summary,
        all_columns=args.all_columns, contig=args.contig, min_base_quality=args.min_base_quality,
        min_depth=args.min_depth, min_minor_frac=args.min_minor_frac,
    )
    print(json.dumps(summary, indent=2))
    return summary
