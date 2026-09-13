"""redwood composite: journal styles, page and label options, legend markup and layout."""
from __future__ import annotations

import random
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest
from matplotlib.text import Text

from redwood.composite import (
    STYLES,
    figure_label,
    panel_letters,
    parse_length_mm,
    parse_markup,
    parse_page,
    plot_journal_figure,
)
from redwood.numts import merge_loci, write_loci


@pytest.fixture
def inputs(tmp_path):
    rng = random.Random(3)
    mito = "".join(rng.choice("ACGT") for _ in range(8000))
    (tmp_path / "mt.fa").write_text(">mt\n" + mito + "\n")
    (tmp_path / "mt.gff").write_text("\n".join("\t".join(f) for f in [
        ["mt", "t", "CDS", "1", "3000", ".", "+", "0", "Name=COX1"],
        ["mt", "t", "tRNA", "3001", "3070", ".", "+", ".", "Name=tRNA-Asn"],
        ["mt", "t", "rRNA", "3071", "4200", ".", "+", ".", "Name=rrnS"],
    ]) + "\n")
    (tmp_path / "nuc.fa").write_text(">chr1\nACGT\n")
    (tmp_path / "nuc.fa.fai").write_text("chr1\t3000000\t6\t4\t5\nchr2\t2000000\t6\t4\t5\n")
    rows = merge_loci([
        {"chrom": "chr1", "tstart": 100000, "tend": 102000, "qstart": 500, "qend": 2500, "strand": "+", "alen": 2000,
         "identity": 0.95, "identity_snv": 0.95, "snv": 100, "indel": 0, "source": "minimap2"},
        {"chrom": "chr2", "tstart": 50000, "tend": 50400, "qstart": 5000, "qend": 5400, "strand": "-", "alen": 400,
         "identity": 0.82, "identity_snv": 0.82, "snv": 72, "indel": 0, "source": "minimap2"},
    ], 8000)
    write_loci(rows, tmp_path / "loci.tsv")
    return tmp_path, rows, {"chr1": 3000000, "chr2": 2000000}


def test_page_and_length_parsing():
    assert parse_page("letter") == (215.9, 279.4) and parse_page("A4") == (210.0, 297.0)
    assert parse_page("210x279mm") == (210.0, 279.0) and parse_page("8.5x11in") == pytest.approx((215.9, 279.4))
    assert parse_length_mm("6.7in") == pytest.approx(170.18) and parse_length_mm("17cm") == 170.0 and parse_length_mm(170) == 170.0
    with pytest.raises(ValueError):
        parse_page("huge")


def test_labels_prefixes_and_panel_letters():
    nc, gen = STYLES["nature-communications"], STYLES["generic"]
    assert figure_label(nc, number=4, supplementary=True) == "Supplementary Fig. 4"
    assert figure_label(nc, number=2) == "Fig. 2"
    assert figure_label(gen, number=4, supplementary=True) == "Supplementary Figure 4"
    assert figure_label(nc, number="S4", prefix="Figure") == "Figure S4"
    assert figure_label(nc) == ""
    assert panel_letters(4) == ["a", "b", "c", "d"] and panel_letters(3, "upper") == ["A", "B", "C"]


def test_markup_keeps_punctuation_against_styled_words():
    units = parse_markup("**Fig. 1 |** genome of *Genus species*. Done")
    assert units[0] == [("Fig.", True, False)] and units[2] == [("|", True, False)]
    assert [u[0][0] for u in units] == ["Fig.", "1", "|", "genome", "of", "Genus", "species", "Done"]
    assert units[5] == [("Genus", False, True)]
    assert units[6] == [("species", False, True), (".", False, False)]
    assert units[-1] == [("Done", False, False)]


def _texts(fig):
    return [t for t in fig.findobj(Text) if t.get_text().strip() and t.get_visible()]


def test_figure_layout_is_the_figure_width_and_keeps_text_in_range(inputs):
    d, rows, lengths = inputs
    res = plot_journal_figure(mito_fasta=d / "mt.fa", loci=rows, chrom_lengths=lengths, out_base=d / "fig", gff=d / "mt.gff",
                              layout="figure", dpi=50, fileforms=("png",), keep_figure=True, species="Genus species")
    fig = res.pop("figure")
    assert fig.get_size_inches()[0] == pytest.approx(170 / 25.4)          # exactly the figure width: places at 100 % in a 170 mm column
    sizes = {round(t.get_fontsize(), 2) for t in _texts(fig)}
    assert min(sizes) >= 5.0 and max(sizes) <= 8.0                     # 5-7 pt text, 8 pt panel letters
    assert {"a", "b", "c", "d"} <= {t.get_text() for t in _texts(fig)}
    md = Path(res["legend_markdown"]).read_text()
    assert md.startswith("**Mitochondrial genome and nuclear mitochondrial insertions (NUMTs) of *Genus species*.**")
    assert "Circular map of the 8,000-bp mitogenome" in md
    plt.close(fig)


def test_page_layout_two_columns_page_size_prefix_and_uppercase(inputs):
    d, rows, lengths = inputs
    nc = STYLES["nature-communications"]
    res = plot_journal_figure(mito_fasta=d / "mt.fa", loci=rows, chrom_lengths=lengths, out_base=d / "page", gff=d / "mt.gff",
                              layout="page", page="a4", panel_labels="upper", label=figure_label(nc, number=4, supplementary=True),
                              species="Genus species", caption_append="Extra *note*.", dpi=50, fileforms=("png",), keep_figure=True)
    fig = res.pop("figure")
    assert tuple(fig.get_size_inches()) == pytest.approx((210 / 25.4, 297 / 25.4))
    texts = _texts(fig)
    assert {"A", "B", "C", "D"} <= {t.get_text() for t in texts} and not ({"a", "b"} & {t.get_text() for t in texts if t.get_fontsize() == 8.0})
    legend = [t for t in texts if t.get_fontsize() == 7.0 and t.get_va() == "baseline"]
    xs = sorted({round(t.get_position()[0] * 210, 0) for t in legend})
    assert any(x < 20 for x in xs) and any(x > 105 for x in xs)        # words in both 90-mm columns (x in mm from the left edge)
    assert res["layout_mm"]["legend_columns"] == 2 and res["layout_mm"]["figure_width"] == 170.0
    assert Path(res["legend_markdown"]).read_text().startswith("**Supplementary Fig. 4 | Mitochondrial genome")
    assert "Extra *note*." in Path(res["legend_markdown"]).read_text()
    plt.close(fig)


def test_cli_composite(inputs):
    from redwood.cli import main

    d, rows, lengths = inputs
    rc = main(["composite", "--mito-fasta", str(d / "mt.fa"), "--numt-loci", str(d / "loci.tsv"), "--nuclear-fasta", str(d / "nuc.fa"),
               "--gff", str(d / "mt.gff"), "--output-base", str(d / "cli" / "S4"), "--figure-label", "Figure S4",
               "--species", "Genus species", "--page", "letter", "--dpi", "40", "--fileform", "pdf", "png"])
    assert rc == 0
    assert (d / "cli" / "S4.pdf").exists() and (d / "cli" / "S4.png").exists()
    assert (d / "cli" / "S4.legend.md").read_text().startswith("**Figure S4 | Mitochondrial genome")
