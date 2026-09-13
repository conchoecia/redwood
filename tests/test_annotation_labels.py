from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from redwood.renderer import (
    PLOT_LIMIT,
    add_feature_label,
    clean_feature_name,
    collapse_locus_features,
    parse_gff,
    trna_short_label,
)


def test_clean_feature_name_strips_annotator_decorations():
    assert clean_feature_name("COX1 CDS 3' Partial CDS") == "COX1"
    assert clean_feature_name("tRNA-Asn gene") == "tRNA-Asn"
    assert clean_feature_name("rrnL rRNA") == "rrnL"
    assert clean_feature_name("ND1 CDS Partial CDS") == "ND1"
    assert clean_feature_name("  ATP8  ") == "ATP8"
    assert clean_feature_name("ND4L-2 gene (second copy of the 3817 bp tandem duplication)") == "ND4L-2"
    assert clean_feature_name("rrnS-frag rRNA (first 196 bp of rrnS)") == "rrnS-frag"


def test_trna_short_label():
    assert trna_short_label("tRNA-Asn") == "N"
    assert trna_short_label("tRNA-Ser2") == "S2"
    assert trna_short_label("tRNA-Leu(UUR)") == "L"
    assert trna_short_label("trnK") == "K"
    assert trna_short_label("trnL2") == "L2"
    assert trna_short_label("weird") == "weird"


def _feature(ftype, start, stop, name, strand="+"):
    return {"type": ftype, "start": start, "stop": stop, "strand": strand, "name": name}


def test_collapse_keeps_one_record_per_locus():
    feats = [
        _feature("gene", 0, 1536, "COX1"),
        _feature("CDS", 0, 1536, "COX1"),
        _feature("gene", 1537, 1599, "tRNA-Asn"),
        _feature("tRNA", 1537, 1599, "tRNA-Asn"),
        _feature("gene", 2000, 3000, "orphan"),          # gene without a specific record: kept
        _feature("gene", 4000, 5000, "ND5"),
        _feature("mRNA", 4000, 5000, "ND5"),
        _feature("exon", 4000, 4400, "ND5"),
        _feature("exon", 4600, 5000, "ND5"),
        _feature("CDS", 4000, 4400, "ND5"),               # spliced CDS: both parts kept, gene dropped
        _feature("CDS", 4600, 5000, "ND5"),
    ]
    kept = collapse_locus_features(feats)
    kinds = sorted((f["type"], f["start"]) for f in kept)
    assert kinds == [("CDS", 0), ("CDS", 4000), ("CDS", 4600), ("gene", 2000), ("tRNA", 1537)]


def test_parse_gff_collapses_mitofinder_style_records(tmp_path):
    gff = tmp_path / "a.gff"
    gff.write_text(
        "\n".join(
            [
                "mt\tmitofinder\tsource\t1\t16000\t.\t+\t.\tName=source",
                "mt\tmitofinder\tgene\t1\t1536\t.\t+\t0\tName=COX1 gene 3' Partial CDS",
                "mt\tmitofinder\tCDS\t1\t1536\t.\t+\t0\tName=COX1 CDS 3' Partial CDS",
                "mt\tmitofinder\tgene\t1538\t1599\t.\t+\t0\tName=tRNA-Asn gene",
                "mt\tmitofinder\ttRNA\t1538\t1599\t.\t+\t0\tName=tRNA-Asn",
                "mt\tmitofinder\tgene\t1600\t2300\t.\t+\t0\tName=rrnS gene",
                "mt\tmitofinder\trRNA\t1600\t2300\t.\t+\t0\tName=rrnS rRNA",
            ]
        )
        + "\n"
    )
    feats = parse_gff(gff)
    assert [(f["type"], f["name"]) for f in feats] == [("CDS", "COX1"), ("tRNA", "tRNA-Asn"), ("rRNA", "rrnS")]


def _texts(ax):
    return [t.get_text() for t in ax.texts]


def test_short_features_and_trnas_get_labels():
    fig, ax = plt.subplots(figsize=(5.8, 5.8), dpi=80)
    ax.set_xlim(-PLOT_LIMIT, PLOT_LIMIT)
    ax.set_ylim(-PLOT_LIMIT, PLOT_LIMIT)
    length = 18685
    # long gene: label along the arc, no leader line
    add_feature_label(ax, _feature("CDS", 0, 1536, "COX1"), length, 1.007, "#fff", outer_radius=1.178, outer_color="#000")
    n_lines_after_long = len(ax.lines)
    assert _texts(ax) == ["COX1"] and n_lines_after_long == 0
    # 159 bp gene: does not fit along the arc -> radial label outside with a leader
    add_feature_label(ax, _feature("CDS", 2338, 2497, "ATP8"), length, 1.007, "#fff", outer_radius=1.178, outer_color="#000")
    assert _texts(ax) == ["COX1", "ATP8"] and len(ax.lines) == 1
    # tRNA one-letter code
    add_feature_label(ax, _feature("tRNA", 1537, 1599, "tRNA-Asn"), length, 1.094, "#fff", outer_radius=1.178,
                      outer_color="#000", label_text=trna_short_label("tRNA-Asn"), fontsize=4.2, min_fontsize=3.4)
    assert "N" in _texts(ax)
    # without an outer radius a non-fitting label is simply skipped (old behavior)
    add_feature_label(ax, _feature("CDS", 2338, 2497, "ATP8"), length, 1.007, "#fff")
    assert _texts(ax).count("ATP8") == 1
    plt.close(fig)
