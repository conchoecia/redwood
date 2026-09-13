import json
from pathlib import Path

import pysam

from redwood.variants import column_variants, flagged_rows, write_variant_table, write_variants

REF = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"  # 60 bp


def _bam(tmp_path: Path, name: str, ref_copies: int, reads: list[tuple[str, int, str, str]]) -> Path:
    """reads: (qname, ref_start, cigar, seq)"""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"LN": len(REF) * ref_copies, "SN": "mt"}]}
    path = tmp_path / f"{name}.bam"
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for qname, start, cigar, seq in sorted(reads, key=lambda r: r[1]):
            a = pysam.AlignedSegment()
            a.query_name, a.reference_id, a.reference_start, a.mapping_quality = qname, 0, start, 60
            a.cigarstring, a.query_sequence = cigar, seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            out.write(a)
    pysam.index(str(path))
    return path


def test_column_variants_reports_mismatch_insertion_and_deletion(tmp_path):
    exact = REF
    mismatch = REF[:10] + "A" + REF[11:]                      # pos 11 (1-based): ref G -> A
    insertion = REF[:20] + "TT" + REF[20:]                     # 2 bp inserted after column 20
    deletion = REF[:30] + REF[32:]                             # columns 31-32 deleted
    bam = _bam(tmp_path, "single", 1, [
        ("exact", 0, "60M", exact),
        ("mm", 0, "60M", mismatch),
        ("ins", 0, "20M2I40M", insertion),
        ("del", 0, "30M2D28M", deletion),
    ])
    rows, summary = column_variants(bam, REF, min_depth=1, min_minor_frac=0.2, min_base_quality=0)
    assert len(rows) == 60
    r11 = rows[10]
    assert (r11["ref"], r11["G"], r11["A"], r11["depth"]) == ("G", 3, 1, 4)
    assert "minor" in r11["events"] and "mismatch" not in r11["events"]
    r20 = rows[19]
    assert r20["ins"] == 1 and r20["ins_seq"] == "TT" and "insertion" in r20["events"]
    r31 = rows[30]
    assert r31["del"] == 1 and r31["depth"] == 4 and "minor" in r31["events"]
    assert summary["insertion_columns"] == 1 and summary["mismatch_columns"] == 0
    assert summary["minor_allele_columns"] >= 3
    # a column where the majority disagrees with the reference is a mismatch
    rows2, summary2 = column_variants(_bam(tmp_path, "maj", 1, [("a", 0, "60M", mismatch), ("b", 0, "60M", mismatch)]),
                                      REF, min_depth=1, min_base_quality=0)
    assert rows2[10]["major"] == "A" and "mismatch" in rows2[10]["events"] and summary2["mismatch_columns"] == 1


def test_reads_without_base_qualities_are_still_counted(tmp_path):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"LN": len(REF), "SN": "mt"}]}
    path = tmp_path / "noqual.bam"
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for i in range(3):
            a = pysam.AlignedSegment()
            a.query_name, a.reference_id, a.reference_start, a.mapping_quality = f"r{i}", 0, 0, 60
            a.cigarstring, a.query_sequence = "60M", REF
            a.query_qualities = pysam.qualitystring_to_array("!" * 60)   # all zero, as in PacBio CLR fastq
            out.write(a)
    pysam.index(str(path))
    rows, summary = column_variants(path, REF, min_depth=1)
    assert rows[0]["depth"] == 3 and summary["deletion_columns"] == 0 and summary["mismatch_columns"] == 0


def test_column_variants_folds_doubled_reference(tmp_path):
    # read mapped entirely inside the second copy of a doubled reference lands on the same columns
    bam = _bam(tmp_path, "doubled", 2, [("copy1", 0, "60M", REF), ("copy2", 60, "60M", REF[:10] + "A" + REF[11:])])
    rows, _ = column_variants(bam, REF, min_depth=1, min_minor_frac=0.2, min_base_quality=0)
    assert rows[10]["depth"] == 2 and rows[10]["A"] == 1 and len(rows) == 60


def test_write_variants_table_and_summary(tmp_path):
    bam = _bam(tmp_path, "w", 1, [("a", 0, "60M", REF), ("b", 0, "60M", REF[:10] + "A" + REF[11:])])
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">mt\n" + REF + "\n")
    out = tmp_path / "v.tsv"
    summary = write_variants(bam, fasta, out, summary_path=tmp_path / "v.json", min_depth=1, min_minor_frac=0.2,
                             min_base_quality=0)
    lines = out.read_text().splitlines()
    assert lines[0].startswith("pos\tref\tdepth") and len(lines) == 2 and lines[1].startswith("11\tG\t2")
    assert json.loads((tmp_path / "v.json").read_text())["rows_written"] == 1
    assert summary["minor_allele_columns"] == 1
    rows, _ = column_variants(bam, REF, min_depth=1, min_minor_frac=0.2, min_base_quality=0)
    assert write_variant_table(rows, tmp_path / "all.tsv", all_columns=True) == 60
    assert len(flagged_rows(rows)) == 1


def test_cli_variants_on_example_dataset(tmp_path):
    from redwood.cli import main

    repo = Path(__file__).resolve().parents[1]
    dataset = repo / "examples" / "datasets" / "human"
    out = tmp_path / "human.variants.tsv"
    main(["variants", "--mito-fasta", str(dataset / "reference.fa"), "--bam", str(dataset / "reads.mapped.bam"),
          "--output", str(out), "--summary", str(tmp_path / "s.json")])
    summary = json.loads((tmp_path / "s.json").read_text())
    assert summary["columns"] == 16569 and summary["columns_with_depth"] > 0
    assert out.read_text().splitlines()[0].split("\t")[:3] == ["pos", "ref", "depth"]
