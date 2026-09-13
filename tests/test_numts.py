import json
import random
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam
import pytest

from redwood.numts import (
    catalog_numts,
    classify_reads,
    draw_numt_catalog,
    draw_numt_landscape,
    draw_numt_sizes,
    plot_composite_figure,
    junction_support,
    merge_loci,
    mito_coverage,
    plot_numt_figures,
    read_loci,
    write_loci,
)

HAVE_MM2 = bool(shutil.which("minimap2"))


def _rand(n, seed):
    rng = random.Random(seed); return "".join(rng.choice("ACGT") for _ in range(n))


def _mutate(seq, rate, seed):
    rng = random.Random(seed); out = []
    for c in seq:
        out.append(rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c)
    return "".join(out)


def _write(path, records):
    with open(path, "w") as fh:
        for name, seq in records: fh.write(f">{name}\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory):
    d = tmp_path_factory.mktemp("numts")
    mito = _rand(16000, 1)
    full = _mutate(mito, 0.03, 2)                     # a full-length NUMT, ~97 % identical
    frag = _mutate(mito[4000:7000], 0.005, 3)         # a 3 kb fragment, ~99.5 % identical
    chr1 = _rand(60000, 4) + full + _rand(40000, 5)
    chr2 = _rand(30000, 6) + frag + _rand(20000, 7)
    chr3 = _rand(50000, 8)
    _write(d / "mito.fa", [("mt", mito)]); _write(d / "nuclear.fa", [("chr1", chr1), ("chr2", chr2), ("chr3", chr3)])
    return d, mito, {"chr1": len(chr1), "chr2": len(chr2), "chr3": len(chr3)}


def test_merge_loci_classes_and_coverage():
    hits = [
        {"chrom": "c", "tstart": 1000, "tend": 9000, "qstart": 0, "qend": 8000, "strand": "+", "alen": 8000, "identity": 0.98, "identity_snv": 0.985, "snv": 100, "indel": 20, "source": "minimap2"},
        {"chrom": "c", "tstart": 9500, "tend": 17500, "qstart": 8000, "qend": 16000, "strand": "+", "alen": 8000, "identity": 0.96, "identity_snv": 0.97, "snv": 200, "indel": 40, "source": "minimap2"},
        {"chrom": "c", "tstart": 90000, "tend": 90400, "qstart": 100, "qend": 500, "strand": "-", "alen": 400, "identity": 0.90, "identity_snv": 0.9, "snv": 40, "indel": 0, "source": "blastn"},
    ]
    rows = merge_loci(hits, 16000, merge=3000)
    assert [r["class"] for r in rows] == ["full-length", "fragment"]
    assert rows[0]["mito_bp"] == 16000 and rows[0]["n_hits"] == 2 and abs(rows[0]["identity"] - 0.97) < 1e-6
    cov = mito_coverage(rows, 16000)
    assert cov[150] == (151, 2, 0.97) and cov[15000][1] == 1   # locus-level identity


@pytest.mark.skipif(not HAVE_MM2, reason="needs minimap2")
def test_catalog_finds_inserted_numts(synthetic, tmp_path):
    d, mito, lengths = synthetic
    rows, summary = catalog_numts(d / "mito.fa", d / "nuclear.fa", tmp_path, threads=2)
    assert summary["n_loci"] == 2
    full = next(r for r in rows if r["class"] == "full-length"); frag = next(r for r in rows if r["class"] != "full-length")
    assert full["chrom"] == "chr1" and 59000 <= full["start"] <= 61000 and 0.95 < full["identity"] < 0.985
    assert frag["chrom"] == "chr2" and frag["mito_bp"] >= 2800 and frag["identity"] > 0.99
    assert (tmp_path / "numts.loci.tsv").exists() and (tmp_path / "numts.mito_coverage.tsv").exists()
    back = read_loci(tmp_path / "numts.loci.tsv")
    assert [r["locus"] for r in back] == [r["locus"] for r in rows] and back[0]["intervals"] == rows[0]["intervals"]


def test_classify_reads_by_segment_origin(tmp_path):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"LN": 100000, "SN": "chr1"}, {"LN": 16000, "SN": "mt"}]}
    loci = [{"locus": "NUMT001", "chrom": "chr1", "start": 60001, "end": 76000, "class": "full-length"}]
    path = tmp_path / "seg.bam"
    def rec(name, ref, start, cigar, seq, supp=False, rev=False):
        a = pysam.AlignedSegment(); a.query_name, a.reference_id, a.reference_start, a.mapping_quality = name, ref, start, 60
        a.cigarstring, a.query_sequence = cigar, seq; a.is_supplementary = supp; a.is_reverse = rev
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq)); return a
    seq = "A" * 3000
    recs = [
        rec("mito", 1, 100, "3000M", seq),
        rec("wrap", 1, 14000, "2000M1000S", seq), rec("wrap", 1, 0, "2000S1000M", seq, supp=True),
        rec("junction", 1, 12000, "2000M1000S", seq), rec("junction", 0, 76000, "2000S1000M", seq, supp=True),
        rec("chimera", 1, 5000, "2000M1000S", seq), rec("chimera", 0, 5000, "2000S1000M", seq, supp=True),
        rec("numt_read", 0, 62000, "3000M", seq),
        rec("nuclear", 0, 10000, "3000M", seq),
        rec("homology", 1, 3000, "2000M1000S", seq), rec("homology", 0, 64000, "2000S1000M", seq, supp=True),   # nuclear part inside the NUMT
        rec("near", 0, 77000, "3000M", seq),                                                                   # nuclear-only, next to the locus
    ]
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for a in sorted(recs, key=lambda a: (a.reference_id, a.reference_start)): out.write(a)
    pysam.index(str(path))
    got = {r["read"]: r["class"] for r in classify_reads(path, "mt", loci)}
    assert got == {"mito": "mito_only", "wrap": "mito_multisegment", "junction": "mito+nuclear_at_NUMT_locus",
                   "chimera": "mito+nuclear_elsewhere", "numt_read": "nuclear_only_at_NUMT_locus", "nuclear": "nuclear_only",
                   "homology": "mito+NUMT_homology", "near": "nuclear_only"}
    js = junction_support(path, loci, flank=500)
    assert js[0]["left_junction_reads"] == 0 and js[0]["right_junction_reads"] == 0   # no read spans a junction by >= 500 bp on both sides here


def test_numt_figures_render(synthetic, tmp_path):
    d, mito, lengths = synthetic
    rows = merge_loci([
        {"chrom": "chr1", "tstart": 60000, "tend": 76000, "qstart": 0, "qend": 16000, "strand": "+", "alen": 16000, "identity": 0.97, "identity_snv": 0.97, "snv": 480, "indel": 0, "source": "minimap2"},
        {"chrom": "chr2", "tstart": 30000, "tend": 33000, "qstart": 4000, "qend": 7000, "strand": "+", "alen": 3000, "identity": 0.995, "identity_snv": 0.995, "snv": 15, "indel": 0, "source": "minimap2"},
    ], 16000)
    outs = plot_numt_figures(rows, 16000, lengths, tmp_path, dpi=50)
    assert [p.name for p in outs] == ["numts.landscape.png", "numts.catalog.png", "numts.sizes.png"]
    assert all(p.exists() and p.stat().st_size > 0 for p in outs)
    fig, (a1, a2, a3) = plt.subplots(3, 1)
    draw_numt_landscape(a1, rows, lengths, min_visible_frac=0.05); draw_numt_catalog(a2, rows, 16000); draw_numt_sizes(a3, rows, 16000)
    assert len(a1.patches) == 3 + 2 and len(a2.lines) >= 2
    # loci are drawn at their true span unless narrower than the visibility minimum (here 5 % of the longest sequence);
    # the 3 kb fragment is widened to that minimum and the title says so, the 16 kb full-length locus keeps its span
    min_w = 0.05 * max(lengths.values())
    widths = sorted(pt.get_width() for pt in a1.patches[3:])
    assert abs(widths[0] - min_w / 1e6) < 1e-9 and abs(widths[1] - 16000 / 1e6) < 1e-9
    assert f"boxes \u2265 {min_w / 1000:.0f} kb wide" in a1.get_title(loc="left") and "wider than real" in a1.get_title(loc="left")
    fig2, a4 = plt.subplots(); draw_numt_landscape(a4, rows, lengths, min_visible_frac=0.0001)
    assert "wider than real" not in a4.get_title(loc="left")   # nothing widened, nothing claimed
    plt.close(fig2)
    # the class legend states the thresholds, so "fragment / large / full-length" is defined on the figure itself
    legend_text = [t.get_text() for t in a2.get_legend().get_texts()]
    assert legend_text == ["full-length (\u226595%)", "large (\u22655 kb)", "fragment (<5 kb)", "alignments in one locus"]
    # one bar per locus on an axis spanning the whole mitogenome, plus the two dashed threshold lines
    bars = [pt for pt in a3.patches if pt.get_width() > 0]
    assert len(bars) == 2 and a3.get_xlim() == (0.0, 16.0)
    assert sorted(round(l.get_xdata()[0], 2) for l in a3.lines if l.get_linestyle() == "--") == [5.0, 15.2]
    plt.close(fig)


def test_composite_figure_is_a_fixed_letter_text_area_page(synthetic, tmp_path):
    from matplotlib.image import imread

    d, mito, lengths = synthetic
    (tmp_path / "mito.fa").write_text(">mt\n" + mito + "\n")
    rows = merge_loci([
        {"chrom": "chr1", "tstart": 60000, "tend": 76000, "qstart": 0, "qend": 16000, "strand": "+", "alen": 16000, "identity": 0.97, "identity_snv": 0.97, "snv": 480, "indel": 0, "source": "minimap2"},
    ], 16000)
    outs = plot_composite_figure(mito_fasta=tmp_path / "mito.fa", loci=rows, chrom_lengths=lengths, out_base=tmp_path / "composite", dpi=40, fileforms=("png",))
    assert outs[0].exists()
    h, w = imread(outs[0]).shape[:2]
    # the page is exactly 6.5 x 9 in (the text area of letter paper with 1-inch margins); nothing may grow or crop it
    assert (w, h) == (260, 360), f"composite is {w / 40:.2f} x {h / 40:.2f} in, expected 6.5 x 9"


@pytest.mark.skipif(not (HAVE_MM2 and shutil.which("samtools")), reason="needs minimap2 and samtools")
def test_cli_numts_end_to_end(synthetic, tmp_path):
    from redwood.cli import main

    d, mito, lengths = synthetic
    # reads: three mito reads, one NUMT-derived read with 1 kb of flank
    nuclear = dict((n, "".join(l.strip() for l in open(d / "nuclear.fa") if not l.startswith(">"))) for n in ["all"])
    chr1 = nuclear["all"][:116000]
    reads = [("m1", mito[1000:9000]), ("m2", mito[5000:15000]), ("m3", mito[12000:16000] + mito[:3000]), ("numt", chr1[59000:70000])]
    with open(tmp_path / "reads.fq", "w") as fq:
        for n, s in reads: fq.write(f"@{n}\n{s}\n+\n{'I' * len(s)}\n")
    out = tmp_path / "numts"
    main(["numts", "--mito-fasta", str(d / "mito.fa"), "--nuclear-fasta", str(d / "nuclear.fa"), "--outdir", str(out),
          "--long-reads", str(tmp_path / "reads.fq"), "--threads", "2", "--dpi", "50", "--figure"])
    summary = json.loads((out / "numts.summary.json").read_text())
    assert summary["n_loci"] == 2 and summary["read_classes"].get("mito_only", 0) >= 2
    assert summary["read_classes"].get("mito+nuclear_at_NUMT_locus", 0) + summary["read_classes"].get("nuclear_only_at_NUMT_locus", 0) >= 1
    assert (out / "numts.read_classes.tsv").exists() and (out / "numts.junction_support.tsv").exists()
    assert (out / "mitogenome_numts.png").exists() and (out / "numts.landscape.png").exists()



def test_catalog_annotation_strip_is_labeled_and_separate(synthetic, tmp_path):
    d, mito, lengths = synthetic
    gff = tmp_path / "mt.gff"
    gff.write_text("\n".join("\t".join(f) for f in [
        ["mt", "t", "CDS", "1", "5000", ".", "+", "0", "Name=COX1"],
        ["mt", "t", "tRNA", "5001", "5070", ".", "+", ".", "Name=tRNA-Asn"],
        ["mt", "t", "rRNA", "5071", "5400", ".", "+", ".", "Name=rrnS"],
        ["mt", "t", "CDS", "9000", "16000", ".", "+", "0", "Name=ND5"],
    ]) + "\n")
    rows = merge_loci([
        {"chrom": "chr2", "tstart": 30000, "tend": 33000, "qstart": 4000, "qend": 7000, "strand": "+", "alen": 3000, "identity": 0.995, "identity_snv": 0.995, "snv": 15, "indel": 0, "source": "minimap2"},
    ], 16000)
    fig, ax = plt.subplots(figsize=(7.5, 3))
    draw_numt_catalog(ax, rows, 16000, gff)
    assert len(ax.child_axes) == 1                                          # the strip is attached to the catalog axes
    strip = ax.child_axes[0]
    assert strip.get_ylabel() == "annotation" and [t.get_text() for t in strip.get_yticklabels()] == ["tRNAs", "genes"]
    names = {t.get_text() for t in strip.texts}
    assert {"COX1", "ND5"} <= names and "non-coding" in names          # big genes named, the 5.4-9 kb gap marked
    assert strip.get_xlabel() == "mitogenome position (kb)" and ax.get_xlabel() == ""
    fig.canvas.draw()
    assert abs(strip.get_position().y1 - ax.get_position().y0) < 1e-3       # strip sits directly under the catalog
    ax.set_position([0.2, 0.5, 0.6, 0.3]); fig.canvas.draw()              # and stays attached when the catalog moves
    assert abs(strip.get_position().y1 - 0.5) < 1e-3 and abs(strip.get_position().x0 - 0.2) < 1e-3
    assert not [pt for pt in ax.patches]                                   # no annotation blocks among the NUMT segments
    plt.close(fig)


def test_repeat_multi_hits_do_not_inflate_the_class():
    # one 500-bp nuclear segment matching 12 copies of a 500-bp mitochondrial repeat unit: 6 kb of mitogenome coverage,
    # but only 500 bp of mtDNA in the locus, so it stays a fragment
    hits = [{"chrom": "chr1", "tstart": 10000, "tend": 10500, "qstart": 2000 + 500 * k, "qend": 2500 + 500 * k, "strand": "+",
             "alen": 500, "identity": 0.97, "identity_snv": 0.97, "snv": 15, "indel": 0, "source": "minimap2"} for k in range(12)]
    (row,) = merge_loci(hits, 16000)
    assert row["mito_bp"] == 6000 and row["nuclear_bp"] == 500 and row["mtdna_bp"] == 500 and row["class"] == "fragment"


def test_read_loci_accepts_tables_without_content_columns(tmp_path):
    old = tmp_path / "old.tsv"
    old.write_text("locus\tchrom\tstart\tend\tspan_bp\tn_hits\tmito_bp\tfrac_of_mito\tidentity\tidentity_snv\tn_snv\tn_indel\tclass\tstrands\tmito_intervals\n"
                   "NUMT001\tchr1\t1\t3000\t3000\t1\t3000\t0.2\t0.99\t0.99\t30\t0\tfragment\t+\t1-3000\n")
    (r,) = read_loci(old)
    assert r["nuclear_bp"] == 3000 and r["mtdna_bp"] == 3000


def _write_bam(path, ref_len, records, ref_name="chr1"):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": ref_name, "LN": ref_len}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, start, cigar, flag in sorted(records, key=lambda r: r[1]):
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.reference_id, a.reference_start, a.cigarstring, a.flag, a.mapping_quality = name, 0, start, cigar, flag, 60
            qlen = sum(n for op, n in a.cigartuples if op in (0, 1, 4, 7, 8))
            a.query_sequence = "A" * qlen; a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
            out.write(a)
    pysam.index(str(path))


def test_junction_support_needs_aligned_bases_across_the_boundary(tmp_path):
    from redwood.numts import junction_support
    bam = tmp_path / "j.bam"
    _write_bam(bam, 100000, [
        ("span", 38000, "4000M", 0),                    # crosses the left boundary into the insertion
        ("span", 39500, "1500M", 2048),                 # supplementary record of the same read: counted once
        ("del", 38000, "2000M3000D2000M", 0),           # a haplotype without the insertion: deletion over the whole locus
    ])
    (row,) = junction_support(bam, [{"locus": "NUMT001", "chrom": "chr1", "start": 40001, "end": 43000, "class": "fragment"}], flank=1000)
    assert row["left_junction_reads"] == 1 and row["right_junction_reads"] == 0
    (strict,) = junction_support(bam, [{"locus": "NUMT001", "chrom": "chr1", "start": 40001, "end": 43000, "class": "fragment"}], flank=1000, min_mapq=61)
    assert strict["left_junction_reads"] == 0                                # all test reads have MAPQ 60


def test_blast_search_fails_loudly_when_blast_is_missing(tmp_path, monkeypatch):
    import redwood.numts as nm
    monkeypatch.setattr(nm.shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError, match="blastn"):
        nm.blastn_hits(tmp_path / "m.fa", tmp_path / "n.fa", tmp_path)


def test_read_loci_of_an_empty_table_is_empty(tmp_path):
    (tmp_path / "empty.tsv").write_text("")
    assert read_loci(tmp_path / "empty.tsv") == []


def test_identity_colorbar_marks_the_saturated_end():
    from redwood.numts import _identity_colorbar
    fig, ax = plt.subplots()
    cb = _identity_colorbar(fig, ax, fraction=0.05, label_size=6, tick_size=5)
    fig.canvas.draw()
    assert cb.extend == "min" and cb.ax.get_yticklabels()[0].get_text() == "\u226475"
    plt.close(fig)
