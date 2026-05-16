"""Workflow helpers for mapping reads and building redwood plots."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pysam

from .renderer import run_plot


MITO_HEADER_TOKENS = (
    "chrmt",
    "chrm",
    " mitochond",
    "mitochondrion",
    "mitochondrial",
    "mt_dna",
    "mtdna",
)


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise SystemExit(f"required tool not found on PATH: {name}")
    return path


def run_command(cmd: list[str], dry_run: bool = False) -> None:
    print(" ".join(cmd), file=sys.stderr)
    if not dry_run:
        subprocess.run(cmd, check=True)


def run_pipeline(cmds: list[list[str]], dry_run: bool = False) -> None:
    print(" | ".join(" ".join(cmd) for cmd in cmds), file=sys.stderr)
    if dry_run:
        return
    processes = []
    previous_stdout = None
    for index, cmd in enumerate(cmds):
        stdin = previous_stdout
        stdout = subprocess.PIPE if index < len(cmds) - 1 else None
        process = subprocess.Popen(cmd, stdin=stdin, stdout=stdout)
        if previous_stdout is not None:
            previous_stdout.close()
        previous_stdout = process.stdout
        processes.append(process)
    return_codes = [process.wait() for process in processes]
    if any(code != 0 for code in return_codes):
        raise subprocess.CalledProcessError(return_codes[-1], " | ".join(" ".join(cmd) for cmd in cmds))


def fasta_records(path: Path):
    name = None
    chunks = []
    opener = open
    if path.suffix == ".gz":
        import gzip

        opener = gzip.open
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def first_fasta_record(path: Path) -> tuple[str, str]:
    try:
        return next(fasta_records(path))
    except StopIteration as exc:
        raise ValueError(f"no FASTA records found in {path}") from exc


def write_fasta(path: Path, name: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start : start + 80] + "\n")


def write_doubled_mito_reference(mito_fasta: Path, output: Path) -> tuple[str, int]:
    name, sequence = first_fasta_record(mito_fasta)
    short_name = name.split()[0]
    write_fasta(output, short_name, sequence + sequence)
    return short_name, len(sequence)


def should_skip_bait_header(header: str, extra_tokens: list[str]) -> bool:
    normalized = header.lower().replace("_", "").replace("-", "")
    first_word = header.split()[0].lower().replace("_", "").replace("-", "")
    if first_word in {"chrm", "chrmt", "mt", "mtdna"}:
        return True
    tokens = list(MITO_HEADER_TOKENS) + [token.lower() for token in extra_tokens]
    return any(token.replace("_", "").replace("-", "") in normalized for token in tokens)


def build_rna_bait_reference(
    mito_fasta: Path,
    nuclear_fasta: Path,
    output: Path,
    extra_exclude_tokens: list[str],
) -> tuple[str, int]:
    mito_name, mito_sequence = first_fasta_record(mito_fasta)
    mito_name = mito_name.split()[0]
    kept = 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as dest:
        for header, sequence in fasta_records(nuclear_fasta):
            if should_skip_bait_header(header, extra_exclude_tokens):
                continue
            kept += 1
            dest.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                dest.write(sequence[start : start + 80] + "\n")
        dest.write(f">{mito_name}\n")
        for start in range(0, len(mito_sequence), 80):
            dest.write(mito_sequence[start : start + 80] + "\n")
    return mito_name, kept


def alignment_depth(read: pysam.AlignedSegment) -> int:
    return sum(stop - start for start, stop in read.get_blocks())


def circular_span(read: pysam.AlignedSegment, sequence_length: int) -> tuple[int, int]:
    start = read.reference_start % sequence_length
    span = read.reference_length or read.query_alignment_length or read.query_length or 0
    return start, min(span, sequence_length)


def add_circular_interval_depth(depth: list[int], start: int, stop: int, sequence_length: int) -> None:
    span = max(0, stop - start)
    if span <= 0:
        return
    if span >= sequence_length:
        for pos in range(sequence_length):
            depth[pos] += 1
        return
    start %= sequence_length
    stop = start + span
    if stop <= sequence_length:
        for pos in range(start, stop):
            depth[pos] += 1
    else:
        for pos in range(start, sequence_length):
            depth[pos] += 1
        for pos in range(0, stop - sequence_length):
            depth[pos] += 1


def select_long_circular_reads(
    input_bam: Path,
    output_bam: Path,
    sequence_length: int,
    target_depth: float,
    min_span_fraction: float,
) -> dict[str, float | int]:
    candidates = []
    with pysam.AlignmentFile(input_bam, "rb") as source:
        header = source.header
        for read in source.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            start, span = circular_span(read, sequence_length)
            if span / sequence_length < min_span_fraction:
                continue
            candidates.append({"read": read, "start": start, "span": span})
    if not candidates:
        raise SystemExit("no mapped reads passed the long-read span threshold")

    target_bases = int(target_depth * sequence_length)
    bins = 96
    coverage = [0] * bins
    selected = []
    selected_bases = 0
    unused = sorted(candidates, key=lambda item: item["span"], reverse=True)
    while unused and selected_bases < target_bases:
        best_index = 0
        best_score = None
        for index, item in enumerate(unused[:2000]):
            start = int(item["start"])
            span = int(item["span"])
            if span >= sequence_length:
                covered_bins = set(range(bins))
            else:
                stride = max(1, sequence_length // bins)
                covered_bins = {
                    int(((start + offset) % sequence_length) / sequence_length * bins)
                    for offset in range(0, span, stride)
                }
            novelty = sum(1 / (1 + coverage[bin_index]) for bin_index in covered_bins)
            score = (novelty, span, -max(coverage[bin_index] for bin_index in covered_bins))
            if best_score is None or score > best_score:
                best_index = index
                best_score = score
        item = unused.pop(best_index)
        selected.append(item["read"])
        selected_bases += int(item["span"])
        stride = max(1, sequence_length // bins)
        for offset in range(0, int(item["span"]), stride):
            coverage[int(((int(item["start"]) + offset) % sequence_length) / sequence_length * bins)] += 1

    output_bam.parent.mkdir(parents=True, exist_ok=True)
    tmp_bam = output_bam.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(tmp_bam, "wb", header=header) as dest:
        for read in selected:
            dest.write(read)
    pysam.sort("-o", str(output_bam), str(tmp_bam))
    pysam.index(str(output_bam))
    tmp_bam.unlink()
    return {
        "input_alignments": len(candidates),
        "selected_alignments": len(selected),
        "selected_mean_depth": round(selected_bases / sequence_length, 2),
        "min_bin_depth": min(coverage),
        "max_bin_depth": max(coverage),
    }


def map_reads(
    reference: Path,
    reads: list[Path],
    output_bam: Path,
    preset: str,
    dry_run: bool = False,
) -> None:
    if not dry_run:
        require_tool("minimap2")
        require_tool("samtools")
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["minimap2", "-a", "-x", preset, str(reference)] + [str(path) for path in reads]
    run_pipeline([cmd, ["samtools", "sort", "-o", str(output_bam), "-"]], dry_run=dry_run)
    if not dry_run:
        run_command(["samtools", "index", str(output_bam)])


def filter_bam_to_reference(input_bam: Path, output_bam: Path, reference_name: str) -> int:
    kept = 0
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(input_bam, "rb") as source:
        header = source.header.to_dict()
        target_sq = [sq for sq in header.get("SQ", []) if sq.get("SN") == reference_name]
        if not target_sq:
            raise SystemExit(f"reference {reference_name!r} was not found in {input_bam}")
        header["SQ"] = target_sq
        with pysam.AlignmentFile(output_bam, "wb", header=header) as dest:
            for read in source.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if source.get_reference_name(read.reference_id) != reference_name:
                    continue
                read.reference_id = 0
                if read.next_reference_id >= 0:
                    mate_name = source.get_reference_name(read.next_reference_id)
                    read.next_reference_id = 0 if mate_name == reference_name else -1
                dest.write(read)
                kept += 1
    pysam.index(str(output_bam))
    return kept


def bam_depth_metrics(bam_path: Path, sequence_length: int) -> dict[str, float | int]:
    depth = [0] * sequence_length
    spans = []
    origin_spanning = 0
    alignments = 0
    with pysam.AlignmentFile(bam_path, "rb") as handle:
        for read in handle.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            alignments += 1
            start, span = circular_span(read, sequence_length)
            spans.append(span)
            if start + span > sequence_length or span >= sequence_length:
                origin_spanning += 1
            for block_start, block_stop in read.get_blocks():
                add_circular_interval_depth(depth, block_start, block_stop, sequence_length)
    covered = sum(1 for value in depth if value > 0)
    mean_depth = round(sum(depth) / sequence_length, 2) if sequence_length else 0
    spans.sort()
    median_span = spans[len(spans) // 2] if spans else 0
    return {
        "alignments": alignments,
        "mean_depth": mean_depth,
        "breadth": round(covered / sequence_length, 4) if sequence_length else 0,
        "min_depth": min(depth) if depth else 0,
        "max_depth": max(depth) if depth else 0,
        "origin_spanning_alignments": origin_spanning,
        "median_span_fraction": round(median_span / sequence_length, 4) if sequence_length else 0,
    }


def write_metrics(args: argparse.Namespace) -> dict[str, object]:
    _, sequence = first_fasta_record(Path(args.mito_fasta))
    metrics: dict[str, object] = {"sequence_length": len(sequence), "tracks": {}}
    if args.long_bam:
        metrics["tracks"]["long_reads"] = bam_depth_metrics(Path(args.long_bam), len(sequence))
    if args.rnaseq_bam:
        metrics["tracks"]["rnaseq"] = bam_depth_metrics(Path(args.rnaseq_bam), len(sequence))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))
    return metrics


def prepare_reference(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    write_doubled_mito_reference(Path(args.mito_fasta), outdir / "mitochondrion.doubled.fa")
    if args.nuclear_fasta:
        build_rna_bait_reference(
            Path(args.mito_fasta),
            Path(args.nuclear_fasta),
            outdir / "rnaseq_bait_plus_mito.fa",
            list(args.exclude_token or []),
        )


def map_long(args: argparse.Namespace) -> dict[str, object]:
    outdir = Path(args.outdir)
    ref = outdir / "references" / "mitochondrion.doubled.fa"
    _, sequence_length = write_doubled_mito_reference(Path(args.mito_fasta), ref)
    raw_bam = outdir / "long_reads.raw.bam"
    output_bam = Path(args.output_bam) if args.output_bam else outdir / "long_reads.redwood.bam"
    map_reads(ref, [Path(path) for path in args.long_reads], raw_bam, args.preset, dry_run=args.dry_run)
    if args.dry_run:
        return {"raw_bam": str(raw_bam), "output_bam": str(output_bam)}
    selection = select_long_circular_reads(
        raw_bam,
        output_bam,
        sequence_length,
        args.target_depth,
        args.min_span_fraction,
    )
    return {"raw_bam": str(raw_bam), "output_bam": str(output_bam), "selection": selection}


def map_rnaseq(args: argparse.Namespace) -> dict[str, object]:
    outdir = Path(args.outdir)
    combined_ref = outdir / "references" / "rnaseq_bait_plus_mito.fa"
    mito_name, bait_contigs = build_rna_bait_reference(
        Path(args.mito_fasta),
        Path(args.nuclear_fasta),
        combined_ref,
        list(args.exclude_token or []),
    )
    raw_bam = outdir / "rnaseq.raw.bam"
    output_bam = Path(args.output_bam) if args.output_bam else outdir / "rnaseq.mito.bam"
    map_reads(combined_ref, [Path(path) for path in args.rnaseq_reads], raw_bam, args.preset, dry_run=args.dry_run)
    if args.dry_run:
        return {"raw_bam": str(raw_bam), "output_bam": str(output_bam), "mito_reference": mito_name}
    kept = filter_bam_to_reference(raw_bam, output_bam, mito_name)
    return {
        "raw_bam": str(raw_bam),
        "output_bam": str(output_bam),
        "mito_reference": mito_name,
        "bait_contigs_kept": bait_contigs,
        "mitochondrial_alignments": kept,
    }


def run_end_to_end(args: argparse.Namespace) -> dict[str, object]:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    results: dict[str, object] = {"outdir": str(outdir)}
    long_bam = None
    rnaseq_bam = None

    if args.long_reads:
        long_args = argparse.Namespace(
            mito_fasta=args.mito_fasta,
            long_reads=args.long_reads,
            outdir=outdir,
            output_bam=outdir / "long_reads.redwood.bam",
            preset=args.long_read_preset,
            target_depth=args.long_read_depth,
            min_span_fraction=args.min_span_fraction,
            dry_run=args.dry_run,
        )
        results["long_reads"] = map_long(long_args)
        long_bam = Path(results["long_reads"]["output_bam"])

    if args.rnaseq_reads:
        if not args.nuclear_fasta:
            raise SystemExit("--nuclear-fasta is required when --rnaseq-reads are provided")
        rna_args = argparse.Namespace(
            mito_fasta=args.mito_fasta,
            nuclear_fasta=args.nuclear_fasta,
            rnaseq_reads=args.rnaseq_reads,
            outdir=outdir,
            output_bam=outdir / "rnaseq.mito.bam",
            preset=args.rnaseq_preset,
            exclude_token=args.exclude_token,
            dry_run=args.dry_run,
        )
        results["rnaseq"] = map_rnaseq(rna_args)
        rnaseq_bam = Path(results["rnaseq"]["output_bam"])

    if not args.dry_run:
        metrics_args = argparse.Namespace(
            mito_fasta=args.mito_fasta,
            long_bam=str(long_bam) if long_bam else None,
            rnaseq_bam=str(rnaseq_bam) if rnaseq_bam else None,
            output=outdir / "redwood.metrics.json",
        )
        results["metrics"] = write_metrics(metrics_args)

        if not args.skip_plot:
            if long_bam is None and args.gff is None:
                raise SystemExit("plotting without long reads requires --gff so the mitochondrial length is known")
            doubled = ["main"] if long_bam else []
            plot_args = argparse.Namespace(
                mito_fasta=str(args.mito_fasta),
                main_bam=str(long_bam) if long_bam else None,
                rnaseq_bam=str(rnaseq_bam) if rnaseq_bam else None,
                gff=str(args.gff) if args.gff else None,
                doubled=doubled,
                dpi=args.dpi,
                fileform=args.fileform,
                interlace=False,
                invert=False,
                log=False,
                no_timestamp=True,
                BASENAME=str(outdir / args.plot_name),
                query=["False"],
                small_start="inside",
                sort="ALNLEN",
                ticks=args.ticks,
                transparent=args.transparent,
                max_reads=args.max_reads,
                dark=False,
                title=None,
                subtitle=None,
                extra_tracks=[],
            )
            run_plot(plot_args)
            results["plot_base"] = str(outdir / args.plot_name)

    manifest = outdir / "redwood.workflow.json"
    manifest.write_text(json.dumps(results, indent=2, default=str) + "\n")
    return results
