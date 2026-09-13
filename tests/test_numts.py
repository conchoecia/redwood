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
    ]
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for a in sorted(recs, key=lambda a: (a.reference_id, a.reference_start)): out.write(a)
    pysam.index(str(path))
    got = {r["read"]: r["class"] for r in classify_reads(path, "mt", loci)}
    assert got == {"mito": "mito_only", "wrap": "mito_multisegment", "junction": "mito+nuclear_at_NUMT_locus",
                   "chimera": "mito+nuclear_elsewhere", "numt_read": "nuclear_only_at_NUMT_locus", "nuclear": "nuclear_only"}
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
    draw_numt_landscape(a1, rows, lengths); draw_numt_catalog(a2, rows, 16000); draw_numt_sizes(a3, rows, 16000)
    assert len(a1.patches) == 3 + 2 and len(a2.lines) >= 2
    # the class legend states the thresholds, so "fragment / large / full-length" is defined on the figure itself
    legend_text = [t.get_text() for t in a2.get_legend().get_texts()]
    assert legend_text == ["full-length (>= 95 % of mtDNA)", "large (>= 5 kb of mtDNA)", "fragment (< 5 kb of mtDNA)"]
    # one bar per locus on an axis spanning the whole mitogenome, plus the two dashed threshold lines
    bars = [pt for pt in a3.patches if pt.get_width() > 0]
    assert len(bars) == 2 and a3.get_xlim() == (0.0, 16.0)
    assert sorted(round(l.get_xdata()[0], 2) for l in a3.lines if l.get_linestyle() == "--") == [5.0, 15.2]
    plt.close(fig)


def test_composite_figure_is_letter_proportioned(synthetic, tmp_path):
    from matplotlib.image import imread

    d, mito, lengths = synthetic
    (tmp_path / "mito.fa").write_text(">mt\n" + mito + "\n")
    rows = merge_loci([
        {"chrom": "chr1", "tstart": 60000, "tend": 76000, "qstart": 0, "qend": 16000, "strand": "+", "alen": 16000, "identity": 0.97, "identity_snv": 0.97, "snv": 480, "indel": 0, "source": "minimap2"},
    ], 16000)
    outs = plot_composite_figure(mito_fasta=tmp_path / "mito.fa", loci=rows, chrom_lengths=lengths, out_base=tmp_path / "composite", dpi=40, fileforms=("png",))
    assert outs[0].exists()
    h, w = imread(outs[0]).shape[:2]
    assert 1.15 <= h / w <= 1.45, f"composite page ratio {h / w:.2f} is not letter-like"


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
