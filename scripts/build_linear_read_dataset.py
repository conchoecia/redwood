#!/usr/bin/env python
"""Build the public Hydra linear example from a bounded real-read FASTQ prefix.

Requires minimap2 and `pip install -e '.[ont]'`. Network access is needed only
for uncached inputs; the README renderer uses the committed fixture offline.
"""
from __future__ import annotations

import argparse
import difflib
import gzip
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from urllib.parse import quote
from urllib.request import urlopen

import pysam

from build_real_read_datasets import fetch_text, parse_fasta
from redwood.ont import reprocess_ont_bam, reverse_complement


ACCESSION = "NC_010214.1"
RUN = "SRR18362010"
FASTQ_URL = "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR183/010/SRR18362010/SRR18362010_1.fastq.gz"


def download_prefix(path, records):
    # Fast compression matters for long WGS reads; only mapped molecules will
    # be committed. The partially written download is never treated as a cache.
    with urlopen(FASTQ_URL, timeout=90) as response, gzip.GzipFile(fileobj=response) as source, \
            gzip.open(path, "wt", compresslevel=1) as output:
        for i in range(records):
            lines = [source.readline().decode("utf-8") for _ in range(4)]
            if not all(lines) or not lines[0].startswith("@") or not lines[2].startswith("+"):
                raise ValueError(f"Incomplete or invalid FASTQ record {i + 1}")
            if len(lines[1].strip()) != len(lines[3].strip()):
                raise ValueError(f"Sequence/quality length mismatch at record {i + 1}")
            output.writelines(lines)
            if (i + 1) % 1000 == 0:
                print(f"Downloaded {i + 1:,} records", flush=True)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def annotation(gff, sequence):
    """Keep RefSeq feature coordinates and add the exact inverted-repeat core."""
    lines = ["##gff-version 3", f"##sequence-region {ACCESSION} 1 {len(sequence)}"]
    for line in gff.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[2] not in {"CDS", "rRNA", "tRNA", "pseudogene"}:
            continue
        attrs = dict(part.split("=", 1) for part in fields[8].split(";") if "=" in part)
        # Keep the partial COX1 copy and its pseudo=true qualifier as a gene
        # arrow; use gene symbols rather than protein accessions as CDS labels.
        if fields[2] == "pseudogene":
            fields[2] = "gene"
        name = attrs.get("product") if fields[2] == "tRNA" else attrs.get("gene")
        attrs["Name"] = name or attrs.get("Name", fields[2])
        attrs.pop("Parent", None)
        fields[8] = ";".join(f"{k}={v}" for k, v in attrs.items())
        lines.append("\t".join(fields))
    window = 2000
    match = difflib.SequenceMatcher(None, sequence[:window],
                                   reverse_complement(sequence[-window:]),
                                   autojunk=False).find_longest_match(0, window, 0, window)
    if match.size < 1000:
        raise ValueError("Expected a long inverted-repeat core in the pinned Hydra reference")
    intervals = [(match.a + 1, match.a + match.size, "+"),
                 (len(sequence) - match.b - match.size + 1, len(sequence) - match.b, "-")]
    note = quote("Exact inverted-repeat core: longest exact match between first 2000 bp and reverse complement of last 2000 bp; not a terminal cap annotation", safe="")
    for i, (start, end, strand) in enumerate(intervals, 1):
        lines.append(f"{ACCESSION}\tredwood-example\trepeat_region\t{start}\t{end}\t.\t{strand}\t.\tID=itr-core-{i};Name=ITR_core;rpt_type=inverted;Note={note}")
    return "\n".join(lines) + "\n", intervals


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("examples/datasets/hydra"))
    parser.add_argument("--cache", type=Path, default=Path("examples/datasets/.cache/hydra"))
    parser.add_argument("--records", type=int, default=10000)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args(argv)
    if args.records < 1 or args.threads < 1:
        parser.error("records and threads must be positive")
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.cache.mkdir(parents=True, exist_ok=True)
    source = {}
    for kind in ("fasta", "gff3", "gbwithparts"):
        path = args.cache / f"{ACCESSION}.{kind}"
        if not path.exists():
            path.write_text(fetch_text(ACCESSION, kind))
        source[kind] = path
    header, sequence = parse_fasta(source["fasta"].read_text())
    if len(sequence) != 16314 or "linear" not in source["gbwithparts"].read_text().splitlines()[0]:
        raise ValueError("Unexpected length/topology for the pinned public reference")
    reference = args.outdir / "reference.fa"
    reference.write_text(f">{ACCESSION}\n" + "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n")
    gff, intervals = annotation(source["gff3"].read_text(), sequence)
    (args.outdir / "annotation.gff").write_text(gff)
    fastq = args.cache / f"{RUN}.first{args.records}.fastq.gz"
    if not fastq.exists():
        partial = fastq.with_suffix(".partial")
        print(f"Downloading the first {args.records:,} records of {RUN}", flush=True)
        download_prefix(partial, args.records)
        partial.replace(fastq)
    record_hash = hashlib.sha256()
    with pysam.FastxFile(str(fastq)) as reads:
        count = 0
        for read in reads:
            record_hash.update(f"@{read.name}\n{read.sequence}\n+\n{read.quality}\n".encode())
            count += 1
        if count != args.records:
            raise ValueError("Cached FASTQ prefix does not contain the requested record count")
    mapper = shutil.which("minimap2")
    if not mapper:
        raise SystemExit("minimap2 must be on PATH")
    cmd = [mapper, "-t", str(args.threads), "-ax", "map-ont", "-Y", "--secondary=no", str(reference), str(fastq)]
    sam = args.cache / "initial.sam"
    print("Mapping the public read prefix", flush=True)
    with sam.open("w") as output, (args.cache / "initial.minimap2.log").open("w") as log:
        subprocess.run(cmd, stdout=output, stderr=log, check=True)
    initial = args.cache / "initial.bam"
    with pysam.AlignmentFile(str(sam), "r") as inp, pysam.AlignmentFile(str(initial), "wb", header=inp.header) as out:
        for read in inp:
            if not read.is_unmapped:
                out.write(read)
    sam.unlink()
    # New work directories preserve previous preprocessing reports on reruns.
    index = 1
    work = args.cache / "ont-processed"
    while work.exists():
        work = args.cache / f"ont-processed-{index}"
        index += 1
    processed = reprocess_ont_bam(initial, reference, work, threads=args.threads)
    bam = args.outdir / "reads.mapped.bam"
    # Retain all mapped primary/supplementary records from the bounded prefix.
    # Evidence uses this full subset; the plot's 30-read limit is separate.
    with pysam.AlignmentFile(processed["output_bam"]) as inp:
        records = [r for r in inp if not r.is_unmapped and not r.is_secondary]
        clean_header = inp.header.to_dict()
    # Tool @PG command lines include local cache paths, which are not needed in
    # the public BAM. Record portable mapping options in the manifest instead.
    clean_header.pop("PG", None)
    with pysam.AlignmentFile(str(bam), "wb", header=clean_header) as out:
        for read in records:
            out.write(read)
    pysam.index(str(bam))
    if not records:
        raise ValueError("No target reads survived mapping")
    (args.outdir / "trimming.tsv").write_text((work / "trimming.tsv").read_text())
    primaries = [r for r in records if not r.is_supplementary]
    aligned = sum(n for r in primaries for op, n in r.cigartuples if op in (0, 7, 8))
    manifest = {
        "name": "hydra", "label": "Hydra oligactis", "topology": "linear",
        "reference_accession": ACCESSION, "fasta_header": header,
        "run_accession": RUN, "bioproject": "PRJNA816482", "platform": "ont",
        "study": "https://doi.org/10.1101/gr.277040.122", "fastq_urls": [FASTQ_URL],
        "reference": reference.name, "annotation": "annotation.gff", "bam": bam.name,
        "bam_index": bam.name + ".bai", "bam_reference": ACCESSION,
        "sequence_length": len(sequence), "doubled_reference": False,
        "records_scanned_per_fastq": [args.records], "fastq_prefix_sha256": sha256(fastq),
        "fastq_prefix_records_sha256": record_hash.hexdigest(),
        "record_checksum_format": "Concatenated @read_name\\nsequence\\n+\\nquality\\n; FASTQ header comments excluded",
        "source_sha256": {k: sha256(v) for k, v in source.items()},
        "fixture_sha256": {p.name: sha256(p) for p in (reference, args.outdir / "annotation.gff", bam)},
        "mapped_reads": len({r.query_name for r in records}),
        "mapped_alignments": len(records), "primary_alignments": len(primaries),
        "primary_aligned_bases": aligned, "mean_primary_depth": round(aligned / len(sequence), 3),
        "mapping_options": "-ax map-ont -Y --secondary=no",
        "minimap2_version": subprocess.check_output([mapper, "--version"], text=True).strip(),
        "selection": "All target-mapped primary and supplementary records from a bounded public FASTQ prefix; display uses 30 terminal-balanced reads",
        "display_terminal_window_bp": 100,
        "ont_preprocessing": {
            "method": "Redwood reprocess_ont_bam with default terminal-window adapters, followed by remapping",
            "audit": "trimming.tsv", "parameters": processed["parameters"],
            "counts": {key: processed[key] for key in (
                "reads", "trimmed_read_count", "removed_bases", "removed_primary_aligned_bases",
                "mapped_after", "unmapped_after", "endpoint_changes")},
            "tools": processed["tools"],
        },
        "repeat_annotation": {
            "method": "Longest exact match between the first 2000 reference bases and reverse complement of the last 2000; exact core only",
            "intervals_1based_inclusive": intervals, "terminal_caps": "Not annotated in this public example",
        },
        "interpretation": "Public reference and read sample are from different source specimens; a plotting fixture, not a new assembly or proof of complete physical termini. No nuclear co-mapping/NUMT classification is performed.",
        "redwood_command": "redwood plot --topology linear --mito-fasta reference.fa --gff annotation.gff --main-bam reads.mapped.bam --linear-layout two-column --max-reads 30 --terminal-window 100 --linear-track depth --linear-track ends --linear-track clips --fileform pdf svg png --no-timestamp -T -o hydra_redwood",
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: manifest[k] for k in ("mapped_reads", "primary_alignments", "mean_primary_depth")}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
