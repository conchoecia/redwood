#!/usr/bin/env python
"""Build compact redwood datasets from real sequencing runs."""

from __future__ import annotations

import argparse
import contextlib
import gzip
import json
import re
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pysam


EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

DATASETS = {
    "beroe": {
        "label": "Beroe forskalii manuscript reads",
        "reference": "MG655622",
        "status": "blocked",
        "bioproject": "PRJNA421807",
        "sra_samples": ["SRS2786396", "SRS2786397", "SRS2786398", "SRS2786399", "SRS5502111"],
        "note": (
            "PMC6991124 points to PRJNA421807, but public run accessions were not "
            "returned by NCBI SRA runinfo/search or ENA run search on 2026-05-15."
        ),
    },
    "human": {
        "label": "Homo sapiens rCRS",
        "reference": "NC_012920.1",
        "run": "DRR165688",
        "platform": "ont",
        "fastq": ["https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR165/DRR165688/DRR165688_1.fastq.gz"],
        "batch_records": 1000,
        "max_records": 8000,
        "target_depth": 75,
        "double_reference": True,
        "min_span_fraction": 0.35,
        "overcollect": 4,
    },
    "drosophila": {
        "label": "Drosophila melanogaster",
        "reference": "NC_024511.2",
        "run": "SRR12187559",
        "platform": "ont",
        "fastq": ["https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR121/059/SRR12187559/SRR12187559_1.fastq.gz"],
        "batch_records": 20000,
        "max_records": 160000,
        "target_depth": 75,
        "double_reference": True,
        "min_span_fraction": 0.25,
        "overcollect": 5,
    },
    "mouse": {
        "label": "Mus musculus",
        "reference": "NC_005089.1",
        "run": "DRR188136",
        "platform": "ont",
        "fastq": ["https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR188/DRR188136/DRR188136_1.fastq.gz"],
        "batch_records": 7500,
        "max_records": 30000,
        "target_depth": 75,
        "double_reference": True,
        "min_span_fraction": 0.35,
        "overcollect": 4,
    },
    "sponge": {
        "label": "Ephydatia muelleri",
        "reference": "NC_010202.1",
        "run": "SRR10983242",
        "platform": "pacbio",
        "fastq": ["https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/042/SRR10983242/SRR10983242_subreads.fastq.gz"],
        "batch_records": 10000,
        "max_records": 40000,
        "target_depth": 75,
        "double_reference": True,
        "min_span_fraction": 0.25,
        "overcollect": 3,
    },
}


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise SystemExit(f"required tool not found on PATH: {name}")
    return path


def fetch_text(accession: str, rettype: str) -> str:
    params = urllib.parse.urlencode(
        {"db": "nuccore", "id": accession, "rettype": rettype, "retmode": "text"}
    )
    with urllib.request.urlopen(f"{EFETCH}?{params}", timeout=60) as handle:
        return handle.read().decode("utf-8")


def parse_fasta(text: str) -> tuple[str, str]:
    header = ""
    seq = []
    for line in text.splitlines():
        if line.startswith(">"):
            header = line[1:].strip()
        elif line.strip():
            seq.append(line.strip())
    if not header:
        raise ValueError("no FASTA header returned")
    return header, "".join(seq).upper()


def write_fasta(path: Path, name: str, seq: str) -> None:
    with path.open("w") as handle:
        handle.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            handle.write(seq[i : i + 80] + "\n")


def attr_value(attrs: str, *keys: str) -> str | None:
    parsed = {}
    for item in attrs.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = urllib.parse.unquote(value)
    for key in keys:
        if parsed.get(key):
            return parsed[key]
    return None


def normalize_gff(accession: str, text: str, seqlen: int) -> str:
    rows = [
        "##gff-version 3",
        f"{accession}\tredwood\tregion\t1\t{seqlen}\t.\t+\t.\tIs_circular=true;Name={accession}",
    ]
    seen = set()
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            continue
        _, source, feat_type, start, stop, _, strand, _, attrs = fields
        if feat_type not in {"CDS", "rRNA", "tRNA"}:
            continue
        name = attr_value(attrs, "standard_name", "gene", "product", "Name", "locus_tag", "ID") or feat_type
        name = re.sub(r"^MT-?", "", name)
        out_type = "gene" if feat_type == "CDS" else feat_type
        key = (out_type, start, stop, strand, name)
        if key in seen:
            continue
        seen.add(key)
        rows.append(
            "\t".join([accession, source or "NCBI", out_type, start, stop, ".", strand, ".", f"Name={name}"])
        )
    return "\n".join(rows) + "\n"


def copy_fastq_prefix(url: str, out: Path, max_records: int) -> int:
    count = 0
    with urllib.request.urlopen(url, timeout=120) as response:
        with gzip.GzipFile(fileobj=response) as gz, gzip.open(out, "wt") as dest:
            while count < max_records:
                record = [gz.readline().decode("utf-8") for _ in range(4)]
                if not record[0]:
                    break
                dest.writelines(record)
                count += 1
    return count


def run_command(cmd: list[str], cwd: Path | None = None) -> None:
    print(" ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, cwd=cwd, check=True)


def mapped_count(bam: Path) -> int:
    result = subprocess.run(["samtools", "view", "-c", str(bam)], check=True, text=True, capture_output=True)
    return int(result.stdout.strip())


def mapped_bases(bam: Path) -> int:
    total = 0
    with pysam.AlignmentFile(bam, "rb") as handle:
        for read in handle.fetch(until_eof=True):
            if not read.is_unmapped:
                total += read.reference_length or read.query_alignment_length or read.query_length or 0
    return total


def circular_span(read: pysam.AlignedSegment, sequence_length: int) -> tuple[int, int]:
    start = read.reference_start % sequence_length
    span = read.reference_length or read.query_alignment_length or read.query_length or 0
    return start, min(span, sequence_length)


def selected_read_metrics(bam: Path, sequence_length: int) -> dict[str, object]:
    spans = []
    bins = [0] * 48
    with pysam.AlignmentFile(bam, "rb") as handle:
        for read in handle.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            start, span = circular_span(read, sequence_length)
            spans.append(span)
            first_bin = int(start / sequence_length * len(bins))
            last_bin = int(((start + span - 1) % sequence_length) / sequence_length * len(bins))
            if span >= sequence_length:
                covered = range(len(bins))
            elif start + span <= sequence_length:
                covered = range(first_bin, last_bin + 1)
            else:
                covered = list(range(first_bin, len(bins))) + list(range(0, last_bin + 1))
            for bin_index in covered:
                bins[bin_index] += 1
    if not spans:
        return {"min_depth_bin": 0, "max_depth_bin": 0, "median_span_fraction": 0, "max_span_fraction": 0}
    spans = sorted(spans)
    return {
        "min_depth_bin": min(bins),
        "max_depth_bin": max(bins),
        "median_span_fraction": round(spans[len(spans) // 2] / sequence_length, 3),
        "max_span_fraction": round(max(spans) / sequence_length, 3),
    }


def merge_bams(inputs: list[Path], out: Path) -> None:
    samtools = require_tool("samtools")
    if len(inputs) == 1:
        shutil.copyfile(inputs[0], out)
    else:
        run_command([samtools, "merge", "-f", str(out), *[str(path) for path in inputs]])
    run_command([samtools, "index", str(out)])


def select_long_circular_reads(
    inputs: list[Path],
    out: Path,
    sequence_length: int,
    target_depth: float,
    min_span_fraction: float,
) -> None:
    candidates = []
    header = None
    for bam in inputs:
        with pysam.AlignmentFile(bam, "rb") as handle:
            if header is None:
                header = handle.header
            for read in handle.fetch(until_eof=True):
                if read.is_unmapped:
                    continue
                start, span = circular_span(read, sequence_length)
                if span / sequence_length < min_span_fraction:
                    continue
                candidates.append(
                    {
                        "read": read,
                        "start": start,
                        "span": span,
                        "span_fraction": span / sequence_length,
                    }
                )
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
            covered_bins = {
                int(((start + offset) % sequence_length) / sequence_length * bins)
                for offset in range(0, span, max(1, sequence_length // bins))
            }
            if span >= sequence_length:
                covered_bins = set(range(bins))
            novelty = sum(1 / (1 + coverage[bin_index]) for bin_index in covered_bins)
            score = (novelty, item["span_fraction"], -max(coverage[bin_index] for bin_index in covered_bins))
            if best_score is None or score > best_score:
                best_index = index
                best_score = score
        item = unused.pop(best_index)
        selected.append(item["read"])
        selected_bases += int(item["span"])
        start = int(item["start"])
        span = int(item["span"])
        for offset in range(0, span, max(1, sequence_length // bins)):
            coverage[int(((start + offset) % sequence_length) / sequence_length * bins)] += 1

    tmp = out.with_suffix(".unsorted.bam")
    assert header is not None
    with pysam.AlignmentFile(tmp, "wb", header=header) as handle:
        for read in selected:
            handle.write(read)
    pysam.sort("-o", str(out), str(tmp))
    pysam.index(str(out))
    tmp.unlink()


def align_fastqs(spec: dict[str, object], ref: Path, fastq_paths: list[Path], bam: Path) -> None:
    minimap2 = require_tool("minimap2")
    bwa = require_tool("bwa")
    samtools = require_tool("samtools")
    if spec["platform"] == "ont":
        cmd = [minimap2, "-a", "-x", "map-ont", str(ref), str(fastq_paths[0])]
    elif spec["platform"] == "pacbio":
        cmd = [minimap2, "-a", "-x", "map-pb", str(ref), str(fastq_paths[0])]
    else:
        if not (ref.with_suffix(ref.suffix + ".bwt")).exists():
            run_command([bwa, "index", str(ref)])
        cmd = [bwa, "mem", str(ref), str(fastq_paths[0]), str(fastq_paths[1])]

    align = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    view = subprocess.Popen([samtools, "view", "-b", "-F", "4", "-"], stdin=align.stdout, stdout=subprocess.PIPE)
    assert align.stdout is not None
    align.stdout.close()
    sort = subprocess.run([samtools, "sort", "-o", str(bam), "-"], stdin=view.stdout, check=True)
    if view.stdout is not None:
        view.stdout.close()
    if align.wait() != 0 or view.wait() != 0 or sort.returncode != 0:
        raise SystemExit(f"mapping failed for {bam}")


def stream_fastq_batches(urls: list[str], cache: Path, name: str, batch_records: int, max_records: int):
    with contextlib.ExitStack() as stack:
        handles = []
        for url in urls:
            response = stack.enter_context(urllib.request.urlopen(url, timeout=120))
            handles.append(stack.enter_context(gzip.GzipFile(fileobj=response)))
        scanned = 0
        chunk_index = 0
        while scanned < max_records:
            batch_size = min(batch_records, max_records - scanned)
            out_paths = [cache / f"{name}_{chunk_index:03d}_{i + 1}.fastq.gz" for i in range(len(urls))]
            written = [0 for _ in urls]
            out_handles = [stack.enter_context(gzip.open(path, "wt")) for path in out_paths]
            try:
                for _ in range(batch_size):
                    records = []
                    for handle in handles:
                        record = [handle.readline().decode("utf-8") for _ in range(4)]
                        records.append(record)
                    if any(not record[0] for record in records):
                        break
                    for i, record in enumerate(records):
                        out_handles[i].writelines(record)
                        written[i] += 1
            finally:
                for out_handle in out_handles:
                    out_handle.close()
            if min(written) == 0:
                break
            scanned += min(written)
            yield out_paths, written, scanned, chunk_index
            chunk_index += 1


def build_mapping(spec: dict[str, object], outdir: Path, cache: Path, name: str) -> dict[str, object]:
    require_tool("samtools")

    accession = str(spec["reference"])
    header, seq = parse_fasta(fetch_text(accession, "fasta"))
    time.sleep(0.35)
    gff = normalize_gff(accession, fetch_text(accession, "gff3"), len(seq))

    ref = outdir / "reference.fa"
    gff_path = outdir / "annotation.gff"
    write_fasta(ref, name, seq)
    gff_path.write_text(gff)
    if spec.get("double_reference"):
        mapping_ref = cache / f"{name}.doubled.fa"
        mapping_name = f"{name}_doubled"
        write_fasta(mapping_ref, mapping_name, seq + seq)
    else:
        mapping_ref = ref
        mapping_name = name

    bam = outdir / "reads.mapped.bam"
    batch_records = int(spec["batch_records"])
    max_records = int(spec["max_records"])
    target_bases = int(float(spec["target_depth"]) * len(seq))
    overcollect_bases = int(target_bases * float(spec.get("overcollect", 3)))
    chunk_bams = []
    total_bases = 0
    total_mapped = 0
    scanned = 0
    for fastq_paths, copied, scanned, chunk_index in stream_fastq_batches(
        [str(url) for url in spec["fastq"]], cache, name, batch_records, max_records
    ):
        chunk_bam = cache / f"{name}_{chunk_index:03d}.mapped.bam"
        align_fastqs(spec, mapping_ref, fastq_paths, chunk_bam)
        chunk_mapped = mapped_count(chunk_bam)
        chunk_bases = mapped_bases(chunk_bam)
        if chunk_mapped > 0:
            chunk_bams.append(chunk_bam)
            total_mapped += chunk_mapped
            total_bases += chunk_bases
        if total_bases >= overcollect_bases:
            break
    if not chunk_bams:
        raise SystemExit(f"no reads from {spec['run']} mapped to {accession}; increase max_records or choose another run")
    select_long_circular_reads(
        chunk_bams,
        bam,
        len(seq),
        float(spec["target_depth"]),
        float(spec.get("min_span_fraction", 0.2)),
    )
    count = mapped_count(bam)
    bases = mapped_bases(bam)
    metrics = selected_read_metrics(bam, len(seq))
    return {
        "reference": ref.name,
        "bam_reference": mapping_name,
        "doubled_reference": bool(spec.get("double_reference")),
        "annotation": gff_path.name,
        "bam": bam.name,
        "bam_index": bam.name + ".bai",
        "fasta_header": header,
        "sequence_length": len(seq),
        "records_scanned_per_fastq": [scanned for _ in spec["fastq"]],
        "mapped_alignments": count,
        "mapped_bases": bases,
        "approx_depth": round(bases / len(seq), 2),
        "selection": {
            "strategy": "greedy_long_read_circular_coverage",
            "target_depth": float(spec["target_depth"]),
            "min_span_fraction": float(spec.get("min_span_fraction", 0.2)),
            **metrics,
        },
    }


def write_blocked_manifest(spec: dict[str, object], outdir: Path, name: str) -> dict[str, object]:
    manifest = {"name": name, **spec}
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("examples/datasets"))
    parser.add_argument("--dataset", choices=sorted(DATASETS), action="append")
    parser.add_argument("--keep-cache", action="store_true")
    args = parser.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    cache = args.outdir / ".cache"
    cache.mkdir(parents=True, exist_ok=True)

    manifests = []
    for name in args.dataset or sorted(DATASETS):
        spec = DATASETS[name]
        outdir = args.outdir / name
        outdir.mkdir(parents=True, exist_ok=True)
        print(f"building {name}", file=sys.stderr)
        if spec.get("status") == "blocked":
            manifests.append(write_blocked_manifest(spec, outdir, name))
            continue
        generated = build_mapping(spec, outdir, cache, name)
        redwood_command = "redwood plot --gff annotation.gff --main-bam reads.mapped.bam --no-timestamp -o plot --query False"
        if spec.get("double_reference"):
            redwood_command = "redwood plot --gff annotation.gff --main-bam reads.mapped.bam --doubled main --no-timestamp -o plot --query False"
        manifest = {
            "name": name,
            "label": spec["label"],
            "reference_accession": spec["reference"],
            "run_accession": spec["run"],
            "platform": spec["platform"],
            "fastq_urls": spec["fastq"],
            **generated,
            "redwood_command": redwood_command,
        }
        (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        manifests.append(manifest)
    (args.outdir / "manifest.json").write_text(json.dumps(manifests, indent=2) + "\n")
    if not args.keep_cache:
        shutil.rmtree(cache)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
