import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam
import pytest

from redwood.multipass import read_marks
from redwood.renderer import (
    PLOT_LIMIT,
    add_feature_label,
    add_position_labels,
    add_read_marks,
    add_spiral_read,
    add_variant_ring,
    plot_file,
)
from redwood.variants import column_variants, write_variant_table

REF = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"  # 60 bp, REF[10] == "G"
REPO = Path(__file__).resolve().parents[1]
HUMAN = REPO / "examples" / "datasets" / "human"


def _bam(tmp_path, name, copies, reads):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"LN": len(REF) * copies, "SN": "mt"}]}
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


def _axes():
    fig, ax = plt.subplots(figsize=(5.8, 5.8), dpi=80)
    ax.set_xlim(-PLOT_LIMIT, PLOT_LIMIT)
    ax.set_ylim(-PLOT_LIMIT, PLOT_LIMIT)
    return fig, ax


def test_read_marks_reports_mismatch_insertion_and_deletion(tmp_path):
    mismatch = REF[:10] + "A" + REF[11:]
    insertion = REF[:20] + "TT" + REF[20:]
    deletion = REF[:30] + REF[32:]
    bam = _bam(tmp_path, "m", 1, [("mm", 0, "60M", mismatch), ("ins", 0, "20M2I40M", insertion), ("del", 0, "30M2D28M", deletion)])
    got = {a.query_name: read_marks(a, REF, len(REF), min_indel=1) for a in pysam.AlignmentFile(bam)}
    assert got["mm"] == [(10, "A")]
    assert got["ins"] == [(19, "I")]          # placed on the reference base preceding the insertion
    assert got["del"] == [(30, "D")]          # first deleted column
    strict = {a.query_name: read_marks(a, REF, len(REF), min_indel=3) for a in pysam.AlignmentFile(bam)}
    assert strict["ins"] == [] and strict["del"] == []


def test_read_marks_fold_onto_the_single_copy(tmp_path):
    bam = _bam(tmp_path, "d", 2, [("copy2", 60, "60M", REF[:10] + "A" + REF[11:])])
    (a,) = list(pysam.AlignmentFile(bam))
    assert read_marks(a, REF, len(REF)) == [(10, "A")]   # offset is relative to the alignment start; base compared to REF[70 % 60]


def test_add_read_marks_and_variant_ring_draw_expected_artists():
    fig, ax = _axes()
    marks = [(10, "A"), (25, "C"), (19, "I"), (30, "D")]
    n_all = add_read_marks(ax, 0, marks, 60, 0.85)
    assert n_all == 4 and len(ax.patches) == 4
    fig, ax = _axes()
    n_shared = add_read_marks(ax, 0, marks, 60, 0.85, mark_filter={10})   # only column 10 is "shared"; indels always
    assert n_shared == 3
    fig, ax = _axes()
    assert add_read_marks(ax, 0, marks, 60, 0.85, mark_filter={10}, mark_indels=False) == 1
    fig, ax = _axes()
    rows = [
        {"pos": 11, "events": "mismatch", "minor_frac": 0.4, "ins_frac": 0.0},
        {"pos": 21, "events": "insertion", "minor_frac": 0.0, "ins_frac": 0.3},
        {"pos": 31, "events": "minor", "minor_frac": 0.2, "ins_frac": 0.0},
        {"pos": 41, "events": "lowdepth", "minor_frac": 0.0, "ins_frac": 0.0},
    ]
    assert add_variant_ring(ax, rows, 60) == 3
    plt.close("all")


def test_spiral_reads_carry_marks():
    fig, ax = _axes()
    before = len(ax.lines)
    add_spiral_read(ax, 0, 1.5, 60, 0.85, 0.01, 0.12, "#333333", marks=[(5, "A"), (70, "T"), (40, "I")], mark_filter={5, 10})   # offset 70 folds to column 10
    assert len(ax.lines) - before == 1 + 3      # the spiral itself + three marker points
    plt.close(fig)


def test_neighboring_outside_labels_are_staggered():
    fig, ax = _axes()
    f1 = {"type": "tRNA", "start": 1000, "stop": 1062, "strand": "+", "name": "tRNA-Phe"}
    f2 = {"type": "tRNA", "start": 1070, "stop": 1132, "strand": "+", "name": "tRNA-Glu"}
    for f, txt in ((f1, "F"), (f2, "E")):
        add_feature_label(ax, f, 18000, 1.094, "#fff", outer_radius=1.178, outer_color="#000", label_text=txt,
                          fontsize=4.6, min_fontsize=3.8, prefer_outside=True)
    radii = {t.get_text(): (t.get_position()[0] ** 2 + t.get_position()[1] ** 2) ** 0.5 for t in ax.texts}
    assert set(radii) == {"F", "E"} and radii["E"] > radii["F"] + 0.01
    plt.close(fig)


def test_last_position_label_is_suppressed_near_the_origin():
    fig, ax = _axes()
    add_position_labels(ax, 15187, "#666")
    assert "15,000 bp" not in {t.get_text() for t in ax.texts}
    fig, ax = _axes()
    add_position_labels(ax, 16569, "#666")
    assert "15,000 bp" in {t.get_text() for t in ax.texts}
    plt.close("all")


def test_plot_with_mismatch_marks_and_variant_table(tmp_path):
    ref = "".join(l.strip() for l in open(HUMAN / "reference.fa") if not l.startswith(">"))
    rows, _ = column_variants(HUMAN / "reads.mapped.bam", ref, min_minor_frac=0.05)
    table = tmp_path / "variants.tsv"
    write_variant_table(rows, table)
    for mode in ("shared",):
        base = tmp_path / f"human_{mode}"
        plot_file(output_base=str(base), fileforms=["png"], dpi=60, reference_fasta=HUMAN / "reference.fa",
                  gff=HUMAN / "annotation.gff", main_bam=HUMAN / "reads.mapped.bam", no_timestamp=True,
                  read_mismatches=mode, variant_table=table)
        assert (tmp_path / f"human_{mode}.png").stat().st_size > 0


@pytest.mark.skipif(not (shutil.which("minimap2") and shutil.which("samtools")), reason="needs minimap2 and samtools")
def test_run_writes_variant_tables_and_metrics(tmp_path):
    from redwood.cli import main

    reads = tmp_path / "reads.fq"
    with open(reads, "w") as fq:
        for a in pysam.AlignmentFile(HUMAN / "reads.mapped.bam"):
            if a.query_sequence:
                fq.write(f"@{a.query_name}\n{a.query_sequence}\n+\n{'I' * len(a.query_sequence)}\n")
    out = tmp_path / "run"
    main(["run", "--mito-fasta", str(HUMAN / "reference.fa"), "--long-reads", str(reads), "--outdir", str(out),
          "--skip-plot", "--long-read-preset", "map-hifi"])
    table = out / "redwood.variants.long_reads.tsv"
    assert table.exists() and table.read_text().splitlines()[0].startswith("pos\tref\tdepth")
    metrics = json.loads((out / "redwood.metrics.json").read_text())
    v = metrics["tracks"]["long_reads"]["variants"]
    assert v["columns"] == 16569 and "effective_minor_threshold" in v
