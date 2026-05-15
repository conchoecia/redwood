#!/usr/bin/env python
"""Build RNA-seq mitochondrial BAM tracks with a nuclear-genome bait reference.

This script intentionally does not download whole nuclear genomes. Provide a
local whole-genome FASTA for each species; mitochondrial-looking contigs are
excluded while building the bait reference, then the committed mitochondrial
reference is appended as the target contig.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import pysam


RNA_DATASETS = {
    "human": {
        "label": "Homo sapiens RNA-seq",
        "run": "DRR001175",
        "strategy": "RNA-Seq",
        "layout": "SINGLE",
        "platform": "ILLUMINA",
        "fastq": ["https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR001/DRR001175/DRR001175.fastq.gz"],
        "max_records": 5000000,
        "minimap2_preset": "sr",
    },
    "mouse": {
        "label": "Mus musculus RNA-seq",
        "run": "DRR001494",
        "strategy": "RNA-Seq",
        "layout": "SINGLE",
        "platform": "ILLUMINA",
        "fastq": ["https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR001/DRR001494/DRR001494.fastq.gz"],
        "max_records": 3000000,
        "minimap2_preset": "sr",
    },
    "drosophila": {
        "label": "Drosophila melanogaster RNA-seq",
        "run": "DRR016419",
        "strategy": "RNA-Seq",
        "layout": "SINGLE",
        "platform": "ILLUMINA",
        "fastq": ["https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR016419/DRR016419.fastq.gz"],
        "max_records": 1000000,
        "target_depth": 75,
        "minimap2_preset": "sr",
    },
    "sponge": {
        "label": "Ephydatia muelleri RNA-seq",
        "run": "SRR14102585",
        "strategy": "RNA-Seq",
        "layout": "PAIRED",
        "platform": "ILLUMINA",
        "fastq": [
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/085/SRR14102585/SRR14102585_1.fastq.gz",
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/085/SRR14102585/SRR14102585_2.fastq.gz",
        ],
        "max_records": 1000000,
        "minimap2_preset": "sr",
    },
}

MITO_HEADER_TOKENS = (
    "chrmt",
    " mitochond",
    "mitochondrion",
    "mitochondrial",
    "mt_dna",
    "mtdna",
)


def require_tool(name: str) -> None:
    if not shutil.which(name):
        raise SystemExit(f"required tool not found on PATH: {name}")


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


def should_skip_bait_header(header: str, extra_tokens: list[str]) -> bool:
    first_word = header.split()[0].lower().replace("_", "").replace("-", "")
    if first_word in {"chrm", "chrmt", "mt", "mtdna"}:
        return True
    normalized = header.lower().replace("_", "").replace("-", "")
    tokens = list(MITO_HEADER_TOKENS) + [token.lower() for token in extra_tokens]
    return any(token.replace("_", "").replace("-", "") in normalized for token in tokens)


def append_filtered_fasta(source: Path, dest, extra_exclude_tokens: list[str]) -> int:
    kept = 0
    write_record = False
    opener = gzip.open if source.suffix == ".gz" else open
    with opener(source, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                write_record = not should_skip_bait_header(line[1:].strip(), extra_exclude_tokens)
                if write_record:
                    kept += 1
                    dest.write(line)
            elif write_record:
                dest.write(line)
    return kept


def first_fasta_name(path: Path) -> str:
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                return line[1:].strip().split()[0]
    raise ValueError(f"no FASTA header found in {path}")


def build_combined_reference(
    bait_fasta: Path,
    mito_fasta: Path,
    output: Path,
    extra_exclude_tokens: list[str],
) -> tuple[str, int]:
    output.parent.mkdir(parents=True, exist_ok=True)
    mito_name = first_fasta_name(mito_fasta)
    with output.open("w") as dest:
        kept = append_filtered_fasta(bait_fasta, dest, extra_exclude_tokens)
        with mito_fasta.open() as mito:
            shutil.copyfileobj(mito, dest)
    return mito_name, kept


def run_pipeline(cmds: list[list[str]]) -> None:
    processes = []
    previous_stdout = None
    for index, cmd in enumerate(cmds):
        stdin = previous_stdout
        stdout = subprocess.PIPE if index < len(cmds) - 1 else None
        print(" ".join(cmd), file=sys.stderr)
        process = subprocess.Popen(cmd, stdin=stdin, stdout=stdout)
        if previous_stdout is not None:
            previous_stdout.close()
        previous_stdout = process.stdout
        processes.append(process)
    failures = [process.wait() for process in processes]
    if any(code != 0 for code in failures):
        raise subprocess.CalledProcessError(failures[-1], " | ".join(" ".join(cmd) for cmd in cmds))


def map_rnaseq(combined_ref: Path, fastqs: list[Path], preset: str, raw_bam: Path) -> None:
    cmd = ["minimap2", "-ax", preset, str(combined_ref)] + [str(path) for path in fastqs]
    run_pipeline([cmd, ["samtools", "sort", "-o", str(raw_bam), "-"]])


def filter_mito_bam(raw_bam: Path, out_bam: Path, mito_name: str) -> int:
    kept = 0
    header = None
    with pysam.AlignmentFile(raw_bam, "rb") as source:
        header = source.header.to_dict()
        mito_sq = [sq for sq in header.get("SQ", []) if sq.get("SN") == mito_name]
        header["SQ"] = mito_sq
        with pysam.AlignmentFile(out_bam, "wb", header=header) as dest:
            for read in source.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if source.get_reference_name(read.reference_id) != mito_name:
                    continue
                read.reference_id = 0
                if read.next_reference_id >= 0:
                    read.next_reference_id = 0 if source.get_reference_name(read.next_reference_id) == mito_name else -1
                dest.write(read)
                kept += 1
    subprocess.run(["samtools", "index", str(out_bam)], check=True)
    return kept


def alignment_depth(read: pysam.AlignedSegment) -> int:
    return sum(stop - start for start, stop in read.get_blocks())


def downsample_bam_to_depth(bam_path: Path, target_depth: float) -> tuple[int, float]:
    tmp_bam = bam_path.with_suffix(".tmp.bam")
    with pysam.AlignmentFile(bam_path, "rb") as source:
        sequence_length = source.lengths[0]
        target_bases = target_depth * sequence_length
        reads = [read for read in source.fetch(until_eof=True) if not read.is_unmapped]
        reads.sort(
            key=lambda read: hashlib.sha256(
                f"{read.query_name}:{read.reference_start}:{read.flag}".encode()
            ).hexdigest()
        )
        kept = []
        kept_bases = 0
        for read in reads:
            kept.append(read)
            kept_bases += alignment_depth(read)
            if kept_bases >= target_bases:
                break
        with pysam.AlignmentFile(tmp_bam, "wb", header=source.header) as dest:
            for read in kept:
                dest.write(read)
    subprocess.run(["samtools", "sort", "-o", str(bam_path), str(tmp_bam)], check=True)
    subprocess.run(["samtools", "index", str(bam_path)], check=True)
    tmp_bam.unlink(missing_ok=True)
    return len(kept), round(kept_bases / sequence_length, 2)


def write_manifest(
    path: Path,
    dataset: str,
    config: dict[str, object],
    kept: int,
    bait_contigs: int,
    observed_depth: float | None = None,
) -> None:
    manifest = {
        "name": dataset,
        "label": config["label"],
        "run_accession": config["run"],
        "strategy": config["strategy"],
        "layout": config["layout"],
        "platform": config["platform"],
        "fastq_urls": config["fastq"],
        "records_scanned_per_fastq": config["max_records"],
        "bait_contigs_kept": bait_contigs,
        "mapped_mito_alignments": kept,
        "bam": "rnaseq.mapped.bam",
        "bam_index": "rnaseq.mapped.bam.bai",
        "selection": {
            "strategy": "map_to_nuclear_bait_plus_mito_then_keep_primary_mito_alignments",
            "minimap2_preset": config["minimap2_preset"],
        },
    }
    if observed_depth is not None:
        manifest["selection"]["post_filter"] = (
            f"deterministic_downsample_to_approximately_{config['target_depth']}x_mean_mitochondrial_depth"
        )
        manifest["selection"]["observed_mean_depth"] = observed_depth
    path.write_text(json.dumps(manifest, indent=2) + "\n")


def build_dataset(
    dataset: str,
    bait_fasta: Path,
    datasets_dir: Path,
    cache_dir: Path,
    max_records: int | None,
    preset: str | None,
    run: str | None,
    fastq_urls: list[str] | None,
) -> None:
    config = dict(RNA_DATASETS[dataset])
    if max_records is not None:
        config["max_records"] = max_records
    if preset is not None:
        config["minimap2_preset"] = preset
    if run is not None:
        config["run"] = run
    if fastq_urls is not None:
        config["fastq"] = fastq_urls
    dataset_dir = datasets_dir / dataset
    work_dir = cache_dir / "rnaseq" / dataset
    work_dir.mkdir(parents=True, exist_ok=True)
    fastqs = []
    for index, url in enumerate(config["fastq"]):
        fastq = work_dir / f"{config['run']}_{index + 1}.fastq.gz"
        if not fastq.exists():
            copy_fastq_prefix(str(url), fastq, int(config["max_records"]))
        fastqs.append(fastq)
    combined_ref = work_dir / "bait_plus_mito.fa"
    mito_name, bait_contigs = build_combined_reference(
        bait_fasta,
        dataset_dir / "reference.fa",
        combined_ref,
        extra_exclude_tokens=[str(config["run"])],
    )
    raw_bam = work_dir / "rnaseq.raw.bam"
    out_bam = dataset_dir / "rnaseq.mapped.bam"
    map_rnaseq(combined_ref, fastqs, str(config["minimap2_preset"]), raw_bam)
    kept = filter_mito_bam(raw_bam, out_bam, mito_name)
    observed_depth = None
    if config.get("target_depth"):
        kept, observed_depth = downsample_bam_to_depth(out_bam, float(config["target_depth"]))
    write_manifest(dataset_dir / "rnaseq_manifest.json", dataset, config, kept, bait_contigs, observed_depth)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", choices=sorted(RNA_DATASETS))
    parser.add_argument("--bait-fasta", type=Path, required=True)
    parser.add_argument("--datasets-dir", type=Path, default=Path("examples/datasets"))
    parser.add_argument("--cache-dir", type=Path, default=Path("examples/datasets/.cache"))
    parser.add_argument("--max-records", type=int)
    parser.add_argument("--preset", default=None, help="Override the minimap2 preset, e.g. sr or splice.")
    parser.add_argument("--run-accession", default=None, help="Override the configured run accession.")
    parser.add_argument(
        "--fastq-url",
        action="append",
        default=None,
        help="Override configured FASTQ URL. Repeat for paired-end runs.",
    )
    args = parser.parse_args(argv)

    for tool in ["minimap2", "samtools"]:
        require_tool(tool)
    build_dataset(
        args.dataset,
        args.bait_fasta,
        args.datasets_dir,
        args.cache_dir,
        args.max_records,
        args.preset,
        args.run_accession,
        args.fastq_url,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
