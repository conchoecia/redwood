"""Auditable adapter trimming and remapping of target-mapped ONT molecules."""

from __future__ import annotations

import copy
import csv
import gzip
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam


# Published ONT ligation motifs, not an inference of the user's sequencing kit.
# https://github.com/rrwick/Porechop/blob/master/porechop/adapters.py
TOP = "AATGTACTTCGTTCAGTTACGTATTGCT"
BOTTOM = "GCAATACGTAACTGAACGAAGT"


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


DEFAULT_5P = (TOP, reverse_complement(BOTTOM))
DEFAULT_3P = (BOTTOM, reverse_complement(TOP))


def ont_options(args):
    return {key: getattr(args, "ont_" + key, default) for key, default in (
        ("adapter_5p", None), ("adapter_3p", None), ("overlap", 12),
        ("error_rate", 0.1), ("end_window", 100), ("threads", 4),
        ("source_reads", None))}


def _reference(path):
    from .workflow import fasta_records

    records = list(fasta_records(Path(path)))
    if len(records) != 1 or not records[0][1]:
        raise ValueError("ONT reprocessing requires one nonempty target FASTA sequence")
    return records[0][0].split()[0], records[0][1]


def _fastq(handle, name, sequence, quality):
    handle.write(f"@{name}\n{sequence}\n+\n{quality}\n")


def _target_names(bam_path, name, length):
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if name not in bam.references or bam.get_reference_length(name) != length:
            raise ValueError("ONT BAM target name/length does not match the mapping reference")
        # Include a molecule even when only its supplementary/secondary record
        # maps to the target. This is candidate extraction, not NUMT filtering.
        return {r.query_name for r in bam.fetch(until_eof=True)
                if not r.is_unmapped and r.reference_name == name}


def extract_reads(bam_path, names, output, source_reads=None):
    """Recover each complete molecule once in its original sequencing orientation."""
    recovered = {}

    def accept(name, sequence, quality):
        if name not in names:
            return
        value = (sequence.upper(), quality)
        if name in recovered and recovered[name] != value:
            raise ValueError(f"conflicting full sequences/qualities for ONT read {name!r}")
        recovered[name] = value

    if source_reads:
        for path in source_reads:
            with pysam.FastxFile(str(path)) as reads:
                for r in reads:
                    if r.name in names:
                        if r.quality is None:
                            raise ValueError("ONT source reads must be FASTQ with base qualities")
                        accept(r.name, r.sequence, r.quality)
    else:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for r in bam.fetch(until_eof=True):
                if r.query_name not in names:
                    continue
                if r.is_paired:
                    raise ValueError("ONT reprocessing expects single-molecule, unpaired reads")
                if (r.query_sequence is None or r.query_qualities is None or
                        any(op == 5 for op, n in (r.cigartuples or []))):
                    continue
                accept(r.query_name, r.get_forward_sequence(),
                       "".join(chr(q + 33) for q in r.get_forward_qualities()))
    missing = names - recovered.keys()
    if missing:
        raise ValueError(f"cannot recover full sequence and qualities for {len(missing)} ONT reads "
                         f"(e.g. {sorted(missing)[0]}); provide --ont-source-reads with the original FASTQ")
    with gzip.open(output, "wt") as handle:
        for name in sorted(recovered):
            _fastq(handle, name, *recovered[name])
    return len(recovered)


def _alignment_summary(bam_path, target):
    """Primary target spans and clips; also retain supplementary-only mapping status."""
    rows = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.reference_name != target:
                continue
            row = rows.setdefault(r.query_name, {"mapped": True, "start_1based": "", "end_1based": "",
                                                "left_clip_bp": "", "right_clip_bp": "",
                                                "primary_aligned_bases": 0})
            if not r.is_secondary and not r.is_supplementary:
                row.update(start_1based=r.reference_start + 1, end_1based=r.reference_end,
                           left_clip_bp=r.cigartuples[0][1] if r.cigartuples[0][0] == 4 else 0,
                           right_clip_bp=r.cigartuples[-1][1] if r.cigartuples[-1][0] == 4 else 0,
                           primary_aligned_bases=sum(n for op, n in r.cigartuples if op in (0, 7, 8)))
    return rows


def _write_tsv(path, rows):
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def aligned_bases_removed(bam_path, target, audit):
    """Count removed M/=/X bases from primary target records in full-read coordinates."""
    rows = {row["read"]: row for row in audit}
    counts = {name: 0 for name in rows}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.reference_name != target:
                continue
            row = rows[r.query_name]
            left, right = row["removed_5p_bp"], row["removed_3p_bp"]
            if r.is_reverse:
                left, right = right, left
            length = row["original_length_bp"]
            if sum(n for op, n in r.cigartuples if op in (0, 1, 4, 5, 7, 8)) != length:
                counts[r.query_name] = None  # Source FASTQ differs from original BAM representation.
                continue
            counts[r.query_name] = 0
            pos = 0
            for op, n in r.cigartuples:
                if op in (0, 7, 8):
                    counts[r.query_name] += max(0, min(pos + n, left) - pos)
                    counts[r.query_name] += max(0, pos + n - max(pos, length - right))
                if op in (0, 1, 4, 5, 7, 8):
                    pos += n
    return counts


def reprocess_ont_bam(input_bam, reference_fasta, outdir, *, copies=1,
                      adapter_5p=None, adapter_3p=None, overlap=12, error_rate=0.1,
                      end_window=100, threads=4, source_reads=None):
    """Extract target-mapped reads, trim end windows with Cutadapt, and remap.

    The destination must be new. Original and trimmed reads, removed tails,
    commands, tool reports and an indexed BAM are retained there. Unmapped
    results remain in the BAM; no read filtering is performed after trimming.
    """
    from .workflow import require_tool, write_fasta

    if overlap < 3 or end_window < overlap or threads < 1 or not 0 <= error_rate < 1:
        raise ValueError("ONT overlap must be >=3, end-window >= overlap, threads >=1, and error-rate in [0,1)")
    adapters = {"5p": list(adapter_5p or DEFAULT_5P), "3p": list(adapter_3p or DEFAULT_3P)}
    for side, sequences in adapters.items():
        adapters[side] = [s.upper() for s in sequences]
        if any(len(s) < overlap or set(s.upper()) - set("ACGT") for s in sequences):
            raise ValueError("ONT adapters must contain only A/C/G/T and be at least the minimum overlap length")
    if importlib.util.find_spec("cutadapt") is None:
        raise SystemExit('ONT trimming requires Cutadapt: install with pip install "redwood[ont]" '
                         '(or pip install -e ".[ont]" in a checkout). Use --no-ont-trim to skip it in mapping workflows.')
    minimap = require_tool("minimap2")
    name, sequence = _reference(reference_fasta)
    if copies < 1:
        raise ValueError("reference copies must be positive")
    sequence *= copies
    names = _target_names(input_bam, name, len(sequence))
    if not names:
        raise ValueError("no ONT reads map to the target reference")
    outdir = Path(outdir).resolve()
    if outdir.exists():
        raise ValueError(f"ONT output directory already exists: {outdir}; choose a new --ont-workdir or --outdir")
    outdir.mkdir(parents=True)
    mapping_reference = outdir / "mapping_reference.fa"
    write_fasta(mapping_reference, name, sequence)
    original = outdir / "reads.original.fastq.gz"
    extract_reads(input_bam, names, original, source_reads)
    commands = []
    windows = {}
    for side in ("5p", "3p"):
        path = outdir / f"ends.{side}.fastq.gz"
        with gzip.open(path, "wt") as output, pysam.FastxFile(str(original)) as reads:
            for r in reads:
                size = min(end_window, len(r.sequence) // 2)
                seq = r.sequence[:size] if side == "5p" else r.sequence[len(r.sequence) - size:]
                qual = r.quality[:size] if side == "5p" else r.quality[len(r.quality) - size:]
                _fastq(output, r.name, seq, qual)
        trimmed = outdir / f"ends.{side}.trimmed.fastq.gz"
        cmd = [sys.executable, "-m", "cutadapt", "-j", str(threads),
               "-O", str(overlap), "-e", str(error_rate),
               "--json", str(outdir / f"{side}.cutadapt.json"),
               "--info-file", str(outdir / f"{side}.cutadapt.info.tsv"),
               "-o", str(trimmed)]
        for index, adapter in enumerate(adapters[side]):
            cmd += ["-g" if side == "5p" else "-a", f"ont_{side}_{index + 1}={adapter}"]
        cmd.append(str(path))
        commands.append(cmd)
        (outdir / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        with (outdir / f"{side}.cutadapt.log").open("w") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
        with pysam.FastxFile(str(trimmed)) as reads:
            windows[side] = {r.name: (r.sequence, r.quality) for r in reads}
        if windows[side].keys() != names:
            raise ValueError("Cutadapt changed the candidate read IDs")

    audit = []
    trimmed_reads = outdir / "reads.trimmed.fastq.gz"
    with pysam.FastxFile(str(original)) as reads, gzip.open(trimmed_reads, "wt") as output, \
            gzip.open(outdir / "removed_tails.fastq.gz", "wt") as removed:
        for r in reads:
            size = min(end_window, len(r.sequence) // 2)
            left = size - len(windows["5p"][r.name][0])
            right = size - len(windows["3p"][r.name][0])
            if (left < 0 or right < 0 or
                    windows["5p"][r.name] != (r.sequence[left:size], r.quality[left:size]) or
                    windows["3p"][r.name] != (r.sequence[len(r.sequence) - size:len(r.sequence) - right],
                                              r.quality[len(r.quality) - size:len(r.quality) - right])):
                raise ValueError("unexpected sequence/quality changes in Cutadapt output")
            # Avoid losing an adapter-only short candidate; retain and flag it.
            status = "trimmed" if left or right else "unchanged"
            if left + right >= len(r.sequence):
                left = right = 0
                status = "retained_empty_trim_candidate"
            stop = len(r.sequence) - right
            _fastq(output, r.name, r.sequence[left:stop], r.quality[left:stop])
            left_seq, right_seq = r.sequence[:left], r.sequence[stop:]
            if left:
                _fastq(removed, r.name + "/removed_5p", left_seq, r.quality[:left])
            if right:
                _fastq(removed, r.name + "/removed_3p", right_seq, r.quality[stop:])
            audit.append(dict(read=r.name, status=status, original_length_bp=len(r.sequence),
                              trimmed_length_bp=stop - left, removed_5p_bp=left, removed_3p_bp=right,
                              removed_5p_sequence=left_seq, removed_3p_sequence=right_seq))
    overlaps = aligned_bases_removed(input_bam, name, audit)
    for row in audit:
        row["removed_primary_aligned_bases"] = overlaps[row["read"]]
    _write_tsv(outdir / "trimming.tsv", audit)
    output_bam = outdir / "reads.trimmed.remapped.bam"
    cmd = [minimap, "-t", str(threads), "-ax", "map-ont", "-Y", "--secondary=no",
           str(mapping_reference), str(trimmed_reads)]
    commands.append(cmd)
    (outdir / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
    with tempfile.TemporaryDirectory(dir=outdir) as temporary:
        sam = Path(temporary) / "remapped.sam"
        with sam.open("w") as output, (outdir / "minimap2.log").open("w") as log:
            subprocess.run(cmd, check=True, stdout=output, stderr=log)
        pysam.sort("-@", str(threads), "-o", str(output_bam), str(sam))
    pysam.index(str(output_bam))
    with pysam.AlignmentFile(str(output_bam), "rb") as bam:
        remapped_names = {r.query_name for r in bam.fetch(until_eof=True)}
    if remapped_names != names:
        raise ValueError("remapping changed the candidate read IDs")
    before = _alignment_summary(input_bam, name)
    after = _alignment_summary(output_bam, name)
    fields = ("start_1based", "end_1based", "left_clip_bp", "right_clip_bp", "primary_aligned_bases")
    comparison = []
    for read in sorted(names):
        row = dict(read=read, mapped_before=read in before, mapped_after=read in after)
        for stage, values in (("before", before), ("after", after)):
            row.update({f"{stage}_{field}": values.get(read, {}).get(field, "") for field in fields})
        comparison.append(row)
    _write_tsv(outdir / "alignment_comparison.tsv", comparison)
    import cutadapt

    result = dict(input_bam=str(Path(input_bam).resolve()), output_bam=str(output_bam),
                  original_reads=str(original), trimmed_reads=str(trimmed_reads),
                  mapping_reference=str(mapping_reference), reference_name=name,
                  reference_length_bp=len(sequence), reference_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                  source_reads=[str(Path(p).resolve()) for p in (source_reads or [])],
                  parameters=dict(adapters=adapters, minimum_overlap=overlap, error_rate=error_rate,
                                  end_window=end_window, threads=threads),
                  tools=dict(cutadapt=cutadapt.__version__, pysam=pysam.__version__,
                             samtools=pysam.__samtools_version__,
                             minimap2=subprocess.check_output([minimap, "--version"], text=True).strip()),
                  reads=len(names), trimmed_read_count=sum(r["status"] == "trimmed" for r in audit),
                  removed_bases=sum(r["removed_5p_bp"] + r["removed_3p_bp"] for r in audit),
                  removed_primary_aligned_bases=sum(n for n in overlaps.values() if n is not None),
                  aligned_removal_unknown_reads=sum(n is None for n in overlaps.values()),
                  mapped_after=len(after), unmapped_after=len(names - after.keys()),
                  primary_aligned_bases_before=sum(r["primary_aligned_bases"] for r in before.values()),
                  primary_aligned_bases_after=sum(r["primary_aligned_bases"] for r in after.values()),
                  endpoint_changes=sum((before.get(r, {}).get("start_1based"), before.get(r, {}).get("end_1based")) !=
                                       (after.get(r, {}).get("start_1based"), after.get(r, {}).get("end_1based")) for r in names),
                  report=str(outdir / "ont_preprocessing.json"),
                  notes=["Candidate IDs come from any alignment to the target, not a mitochondrial-origin classifier.",
                         "Read sequences and removed tails use original sequencing orientation; alignment coordinates are 1-based inclusive.",
                         "Only end-window adapter matches are trimmed; no quality, homopolymer, or length filtering.",
                         "Matches remove the adapter and sequence toward the outside read end. Internal adapter matches beyond the window are not processed.",
                         "Reference-only remapping does not distinguish mitochondrial reads from NUMTs."])
    Path(result["report"]).write_text(json.dumps(result, indent=2) + "\n")
    print(f"ONT: trimmed {result['trimmed_read_count']:,}/{len(names):,} mapped candidates; "
          f"{result['mapped_after']:,} remapped. Report: {result['report']}", file=sys.stderr)
    return result


def reprocess_existing(args):
    name, sequence = _reference(args.mito_fasta)
    with pysam.AlignmentFile(str(args.main_bam), "rb") as bam:
        if name not in bam.references:
            raise ValueError("ONT BAM does not contain the FASTA target name")
        length = bam.get_reference_length(name)
    if length % len(sequence) or (getattr(args, "topology", "circular") == "linear" and length != len(sequence)):
        raise ValueError("ONT BAM reference length is incompatible with the FASTA/topology")
    return reprocess_ont_bam(args.main_bam, args.mito_fasta, args.outdir,
                            copies=length // len(sequence), **ont_options(args))


def preprocess_plot(args):
    if not args.main_bam or not args.mito_fasta:
        raise ValueError("--reprocess-ont requires --main-bam and --mito-fasta")
    if getattr(args, "topology", "circular") == "linear" and getattr(args, "doubled", []):
        raise ValueError("linear ONT preprocessing cannot use --doubled")
    prepared = copy.copy(args)
    prepared.outdir = getattr(args, "ont_workdir", None) or Path((args.BASENAME or "redwood") + ".ont-preprocess")
    result = reprocess_existing(prepared)
    prepared.main_bam = result["output_bam"]
    prepared.reprocess_ont = False
    prepared.ont_preprocessing = result
    # The renderer needs to know if the retained circular mapping target is tandem.
    if getattr(args, "topology", "circular") == "circular":
        _, sequence = _reference(args.mito_fasta)
        prepared.doubled = [d for d in getattr(args, "doubled", []) if d != "main"]
        if result["reference_length_bp"] > len(sequence):
            prepared.doubled.append("main")
    return prepared
