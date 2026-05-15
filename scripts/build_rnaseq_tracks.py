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
        "max_records": 200000,
        "minimap2_preset": "sr",
    },
    "mouse": {
        "label": "Mus musculus RNA-seq",
        "run": "DRR001494",
        "strategy": "RNA-Seq",
        "layout": "SINGLE",
        "platform": "ILLUMINA",
        "fastq": ["https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR001/DRR001494/DRR001494.fastq.gz"],
        "max_records": 200000,
        "minimap2_preset": "sr",
    },
    "drosophila": {
        "label": "Drosophila melanogaster RNA-seq",
        "run": "DRR016419",
        "strategy": "RNA-Seq",
        "layout": "SINGLE",
        "platform": "ILLUMINA",
        "fastq": ["https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR016419/DRR016419.fastq.gz"],
        "max_records": 200000,
        "minimap2_preset": "sr",
    },
    "sponge": {
        "label": "Ephydatia muelleri RNA-seq",
        "run": "SRR3168560",
        "strategy": "RNA-Seq",
        "layout": "SINGLE",
        "platform": "ILLUMINA",
        "fastq": ["https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR316/000/SRR3168560/SRR3168560.fastq.gz"],
        "max_records": 200000,
        "minimap2_preset": "sr",
    },
}

MITO_HEADER_TOKENS = (
    "chrm",
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
    normalized = header.lower().replace("_", "").replace("-", "")
    tokens = list(MITO_HEADER_TOKENS) + [token.lower() for token in extra_tokens]
    return any(token.replace("_", "").replace("-", "") in normalized for token in tokens)


def append_filtered_fasta(source: Path, dest, extra_exclude_tokens: list[str]) -> int:
    kept = 0
    write_record = False
    with source.open() as handle:
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
        kept = append_filtered_fasta(bait_fasta, dest, extra_exclude_tokens + [mito_name])
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


def write_manifest(path: Path, dataset: str, config: dict[str, object], kept: int, bait_contigs: int) -> None:
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
    path.write_text(json.dumps(manifest, indent=2) + "\n")


def build_dataset(
    dataset: str,
    bait_fasta: Path,
    datasets_dir: Path,
    cache_dir: Path,
    max_records: int | None,
    preset: str | None,
) -> None:
    config = dict(RNA_DATASETS[dataset])
    if max_records is not None:
        config["max_records"] = max_records
    if preset is not None:
        config["minimap2_preset"] = preset
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
    write_manifest(dataset_dir / "rnaseq_manifest.json", dataset, config, kept, bait_contigs)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", choices=sorted(RNA_DATASETS))
    parser.add_argument("--bait-fasta", type=Path, required=True)
    parser.add_argument("--datasets-dir", type=Path, default=Path("examples/datasets"))
    parser.add_argument("--cache-dir", type=Path, default=Path("examples/datasets/.cache"))
    parser.add_argument("--max-records", type=int)
    parser.add_argument("--preset", default=None, help="Override the minimap2 preset, e.g. sr or splice.")
    args = parser.parse_args(argv)

    for tool in ["minimap2", "samtools"]:
        require_tool(tool)
    build_dataset(args.dataset, args.bait_fasta, args.datasets_dir, args.cache_dir, args.max_records, args.preset)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
