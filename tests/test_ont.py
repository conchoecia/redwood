import csv
import importlib.util
from pathlib import Path
import random
import shutil

import pysam
import pytest

from redwood.cli import build_parser
from redwood.ont import (TOP, DEFAULT_3P, _target_names, extract_reads,
                         reprocess_existing, reverse_complement, aligned_bases_removed)
from redwood.workflow import map_long


def write_bam(path, records, length=2000):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "mito", "LN": length}, {"SN": "nuclear", "LN": 5000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, seq, cigar, flag, target in records:
            read = pysam.AlignedSegment(bam.header)
            read.query_name = name
            read.query_sequence = seq
            read.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            read.flag = flag
            read.reference_id = target
            read.reference_start = 0 if target >= 0 else -1
            read.mapping_quality = 60 if target >= 0 else 0
            read.cigarstring = cigar
            bam.write(read)


def read_fastq(path):
    with pysam.FastxFile(str(path)) as reads:
        return {r.name: (r.sequence, r.quality) for r in reads}


def test_extract_complete_molecules_and_original_orientation(tmp_path):
    bam = tmp_path / "reads.bam"
    seq = "ACGTC" * 20
    write_bam(bam, [("reverse", reverse_complement(seq), "5S90M5S", 16, 0),
                    ("reverse", reverse_complement(seq), "20S80M", 16 | 2048, 0),
                    ("cross", seq, "100M", 0, 1),
                    ("cross", seq[:50], "50M50H", 2048, 0),
                    ("nuclear_only", seq, "100M", 0, 1),
                    ("unmapped", seq, None, 4, -1)])
    names = _target_names(bam, "mito", 2000)
    assert names == {"reverse", "cross"}
    output = tmp_path / "original.fq.gz"
    assert extract_reads(bam, names, output) == 2
    assert read_fastq(output) == {name: (seq, "I" * 100) for name in names}


def test_hard_clipped_only_requires_original_fastq(tmp_path):
    bam = tmp_path / "hard.bam"
    seq = "ACGT" * 25
    write_bam(bam, [("hard", seq[:50], "50M50H", 2048, 0)])
    output = tmp_path / "reads.fq.gz"
    with pytest.raises(ValueError, match="ont-source-reads"):
        extract_reads(bam, {"hard"}, output)
    fastq = tmp_path / "full.fastq"
    fastq.write_text("@hard\n" + seq + "\n+\n" + "I" * 100 + "\n")
    extract_reads(bam, {"hard"}, output, [fastq])
    assert read_fastq(output)["hard"][0] == seq


def test_conflicting_read_names_fail(tmp_path):
    bam = tmp_path / "collision.bam"
    write_bam(bam, [("same", "A" * 100, "100M", 0, 0), ("same", "C" * 100, "100M", 0, 0)])
    with pytest.raises(ValueError, match="conflicting"):
        extract_reads(bam, {"same"}, tmp_path / "reads.fq.gz")


def test_removed_aligned_bases_respect_strand_and_hard_clips(tmp_path):
    bam = tmp_path / "reads.bam"
    write_bam(bam, [("reverse", "A" * 100, "5S90M5S", 16, 0),
                    ("hard", "A" * 90, "10H90M", 0, 0)])
    audit = [dict(read="reverse", original_length_bp=100, removed_5p_bp=9, removed_3p_bp=2),
             dict(read="hard", original_length_bp=100, removed_5p_bp=12, removed_3p_bp=0)]
    assert aligned_bases_removed(bam, "mito", audit) == {"reverse": 4, "hard": 2}


requires_tools = pytest.mark.skipif(
    not importlib.util.find_spec("cutadapt") or not shutil.which("minimap2") or not shutil.which("samtools"),
    reason="ONT integration requires cutadapt, minimap2 and samtools")


@requires_tools
def test_real_ont_workflow_preserves_caps_ids_and_internal_sequence(tmp_path):
    rng = random.Random(712)
    genome = "GGGGGGG" + "".join(rng.choices("ACGT", k=3986)) + "CCCCCCC"
    fasta = tmp_path / "mito.fa"
    fasta.write_text(">mito\n" + genome + "\n")
    adapter3 = DEFAULT_3P[1]
    damaged = TOP[:10] + ("A" if TOP[10] != "A" else "C") + TOP[11:]
    molecules = {
        "both": TOP + genome + adapter3,
        "partial": TOP[-18:] + genome + adapter3[:16],
        "error": damaged + genome,
        "caps": "GGG" + genome + "CCCC",
        "untouched": genome,
        "internal": genome[:2000] + TOP + genome[2000:],
        "reverse": reverse_complement(TOP + genome + adapter3),
        "unmapped": "A" * 300,
    }
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("".join(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n" for name, seq in molecules.items()))
    args = build_parser().parse_args(["advanced", "map-long", "--topology", "linear",
        "--mito-fasta", str(fasta), "--long-reads", str(fastq), "--outdir", str(tmp_path / "work"),
        "--ont-threads", "1"])
    result = map_long(args)
    report = result["ont_preprocessing"]
    original, trimmed = read_fastq(report["original_reads"]), read_fastq(report["trimmed_reads"])
    assert set(trimmed) == set(molecules) - {"unmapped"}
    assert original == {name: (seq, "I" * len(seq)) for name, seq in molecules.items() if name != "unmapped"}
    for name in ("both", "partial", "error"):
        assert trimmed[name] == (genome, "I" * len(genome))
    assert trimmed["reverse"][0] == reverse_complement(genome)
    for name in ("caps", "untouched", "internal"):
        assert trimmed[name] == original[name]
    assert report["reads"] == 7 and report["trimmed_read_count"] == 4
    assert report["mapped_after"] == 7 and report["unmapped_after"] == 0
    with pysam.AlignmentFile(result["output_bam"], "rb") as bam:
        assert bam.get_reference_length("mito") == len(genome)
        assert bam.has_index()
    assert fasta.read_text() == ">mito\n" + genome + "\n"
    assert Path(report["report"]).exists()
    with (Path(report["report"]).parent / "trimming.tsv").open() as f:
        rows = {r["read"]: r for r in csv.DictReader(f, delimiter="\t")}
    assert rows["both"]["removed_5p_sequence"] == TOP
    assert rows["both"]["removed_3p_sequence"] == adapter3
    # Reject reruns before modifying the initial mapping or previous products.
    before = Path(result["raw_bam"]).read_bytes()
    with pytest.raises(ValueError, match="already exists"):
        map_long(args)
    assert Path(result["raw_bam"]).read_bytes() == before


@requires_tools
def test_existing_circular_bam_keeps_tandem_reference(tmp_path):
    rng = random.Random(47)
    genome = "".join(rng.choices("ACGT", k=2000))
    fasta = tmp_path / "mito.fa"
    fasta.write_text(">mito\n" + genome + "\n")
    bam = tmp_path / "tandem.bam"
    seq = TOP + genome * 2
    write_bam(bam, [("two_pass", seq, f"{len(TOP)}S4000M", 0, 0),
                    ("false_candidate", TOP + "A" * 80, f"{len(TOP)}S80M", 0, 0)], length=4000)
    args = build_parser().parse_args(["advanced", "reprocess-ont", "--mito-fasta", str(fasta),
        "--main-bam", str(bam), "--outdir", str(tmp_path / "processed"), "--ont-threads", "1"])
    result = reprocess_existing(args)
    assert result["reference_length_bp"] == 4000
    assert read_fastq(result["trimmed_reads"])["two_pass"][0] == genome * 2
    assert result["reads"] == 2 and result["unmapped_after"] == 1
    with pysam.AlignmentFile(result["output_bam"], "rb") as mapped:
        assert any(r.query_name == "false_candidate" and r.is_unmapped for r in mapped.fetch(until_eof=True))
    args.topology = "linear"
    args.outdir = tmp_path / "invalid"
    with pytest.raises(ValueError, match="incompatible"):
        reprocess_existing(args)
    assert not args.outdir.exists()


def test_plot_explicit_ont_preprocessing_reaches_renderer(tmp_path, monkeypatch):
    from redwood.renderer import run_plot

    args = build_parser().parse_args(["plot", "--topology", "linear", "--reprocess-ont",
        "--mito-fasta", "mito.fa", "--main-bam", "original.bam", "-o", str(tmp_path / "figure")])
    result = {"output_bam": "trimmed.bam", "report": "report.json"}
    monkeypatch.setattr("redwood.ont.reprocess_existing", lambda a: result)
    seen = []
    monkeypatch.setattr("redwood.linear.run_linear_plot", lambda a: seen.append(a))
    run_plot(args)
    assert seen[0].main_bam == "trimmed.bam"
    assert seen[0].ont_preprocessing is result
    assert args.main_bam.endswith("original.bam")


@pytest.mark.parametrize("extra", [["--no-ont-trim"], ["--preset", "map-hifi"]])
def test_mapping_opt_out_and_hifi_skip_preprocessing(tmp_path, monkeypatch, extra):
    fasta = tmp_path / "mito.fa"
    fasta.write_text(">mito\n" + "ACGT" * 500 + "\n")
    def unexpected(*args):
        pytest.fail("ONT preprocessing was called")
    monkeypatch.setattr("redwood.workflow._reprocess_mapped_ont", unexpected)
    args = build_parser().parse_args(["advanced", "map-long", "--topology", "linear",
        "--mito-fasta", str(fasta), "--long-reads", "unused.fq", "--outdir", str(tmp_path / "work"),
        "--dry-run", *extra])
    assert map_long(args)["ont_preprocessing"] is None
