"""Regression tests for the review findings: AT windows at the origin, upright labels, terminal insertion marks, insertion
fractions with low-quality bases, separate gene/tRNA key layers and the key font reservation."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam
import pytest

from redwood.multipass import read_marks
from redwood.renderer import PLOT_LIMIT, add_feature_label, at_profile, draw_circular_plot, track_legend_layers
from redwood.variants import column_variants


def _segment(cigar, start=0, qual=40, name="r"):
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "mt", "LN": 300}]})
    a = pysam.AlignedSegment(header)
    a.query_name, a.reference_id, a.reference_start, a.cigarstring, a.mapping_quality = name, 0, start, cigar, 60
    qlen = sum(n for op, n in a.cigartuples if op in (0, 1, 4, 7, 8))
    a.query_sequence = "A" * qlen
    a.query_qualities = pysam.qualitystring_to_array(chr(qual + 33) * qlen)
    return a, header


def test_at_windows_wrap_around_the_origin():
    ref = "A" * 100 + "C" * 100
    prof = at_profile(ref, window=21)
    assert prof[199] == pytest.approx(10 / 21) and prof[0] == pytest.approx(11 / 21) and prof[50] == 1.0


def _label_axes():
    fig, ax = plt.subplots(figsize=(5.8, 5.8))
    ax.set_aspect("equal"); ax.set_xlim(-PLOT_LIMIT, PLOT_LIMIT); ax.set_ylim(-PLOT_LIMIT, PLOT_LIMIT)
    return fig, ax


def _upright(rotation):
    r = rotation % 360
    return r <= 90 or r >= 270


def test_along_arc_labels_on_the_bottom_half_are_upright():
    fig, ax = _label_axes()
    add_feature_label(ax, {"start": 450, "stop": 550, "name": "ND5", "type": "CDS"}, 1000, 1.03, "white")
    assert ax.texts and _upright(ax.texts[-1].get_rotation())
    plt.close(fig)


def test_radial_labels_on_the_left_half_are_upright():
    fig, ax = _label_axes()
    add_feature_label(ax, {"start": 745, "stop": 755, "name": "V", "type": "tRNA"}, 1000, 1.094, "white",
                      outer_radius=1.18, outer_color="black", prefer_outside=True)
    text = ax.texts[-1]
    assert _upright(text.get_rotation()) and text.get_ha() == "right"
    plt.close(fig)


def test_terminal_insertion_is_marked_but_a_soft_clip_is_not():
    ref = "A" * 300
    ins, _ = _segment("10M3I")
    clip, _ = _segment("10M3S")
    assert (9, "I") in read_marks(ins, ref, 300, min_indel=1)
    assert not [m for m in read_marks(clip, ref, 300, min_indel=1) if m[1] == "I"]


def test_insertion_fraction_uses_read_depth_not_quality_filtered_depth(tmp_path):
    ref = "A" * 300
    bam = tmp_path / "v.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "mt", "LN": 300}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i in range(10):
            a, _ = _segment("50M5I50M", start=100, qual=5, name=f"low{i}")
            out.write(pysam.AlignedSegment.from_dict(a.to_dict(), out.header))
        b, _ = _segment("100M", start=100, qual=30, name="high")
        out.write(pysam.AlignedSegment.from_dict(b.to_dict(), out.header))
    pysam.index(str(bam))
    rows, _ = column_variants(bam, ref, min_base_quality=20, noise_multiplier=0)
    row = rows[149]                                    # insertion after 0-based column 149
    assert row["ins"] == 10 and 0 < row["ins_frac"] <= 1 and row["ins_frac"] == pytest.approx(10 / 11, abs=1e-3)


def test_key_lists_genes_and_trnas_separately():
    keys = [k for k, _ in track_legend_layers(has_rnaseq=False, has_at=True, has_variants=False, has_numts=False,
                                              has_multipass=False, has_regular=False, has_genes=False, has_trna=True)]
    assert keys == ["at", "trna"]


def test_key_font_never_exceeds_the_reserved_label_room():
    fig, ax = plt.subplots(figsize=(3, 3))
    draw_circular_plot(ax, length=2000, reference="ACGT" * 500, track_legend=True)
    key_texts = [t for t in ax.texts if t.get_position()[0] > 1.0 and t.get_position()[1] < 0]
    assert key_texts and all(t.get_fontsize() <= ax._redwood_key_font + 1e-9 for t in key_texts)
    fig.canvas.draw()
    right = ax.transData.inverted().transform((max(t.get_window_extent().x1 for t in key_texts), 0))[0]
    assert right <= ax.get_xlim()[1] + 1e-6                # label text stays inside the axes
    plt.close(fig)
