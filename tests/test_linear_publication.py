import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
import pytest

from redwood.cli import build_parser
from redwood.linear import save_linear_figure
from redwood.linear_evidence import Reference, Segment
from redwood.linear_publication import draw_publication_linear, resolve_linear_layout
from redwood.linear_redwood import add_clip_stack, add_gradient_read, read_outline


@pytest.mark.parametrize("journal,widths", [("nature", (89, 183)), ("nature-communications", (88, 180))])
def test_layout_defaults_and_explicit_overrides(journal, widths):
    for layout, width, count in zip(("one-column", "two-column"), widths, (30, 80)):
        args = build_parser().parse_args(["plot", "--topology", "linear", "--linear-layout", layout,
                                         "--publication-journal", journal])
        assert resolve_linear_layout(args)
        assert args.width * 25.4 == pytest.approx(width)
        assert args.max_reads == count
        args = build_parser().parse_args(["plot", "--linear-layout", layout, "--width", "4", "--max-reads", "0"])
        resolve_linear_layout(args)
        assert args.width == 4 and args.max_reads == 0


@pytest.mark.parametrize("layout", ["one-column", "two-column"])
def test_all_annotation_names_are_readable_without_collisions(layout):
    reference = Reference("mito", "ACGT" * 4092)
    features = [dict(name=name, start=start, stop=stop, type=kind, strand=strand, attributes={})
                for name, start, stop, kind, strand in [
                    ("COX1_C", 0, 820, "CDS", "-"), ("tRNA-Met", 1200, 1260, "tRNA", "+"),
                    ("tRNA-Trp", 2300, 2360, "tRNA", "+"), ("rrnL", 1250, 2350, "rRNA", "+"),
                    ("COX2", 2350, 2990, "CDS", "+"), ("ATP8", 2990, 3150, "CDS", "+"),
                    ("ATP6", 3140, 3740, "CDS", "+"), ("COX3", 3740, 4490, "CDS", "+"),
                    ("ND2", 4490, 5540, "CDS", "+"), ("ND5", 5540, 8000, "CDS", "+"),
                    ("rrnS", 8000, 8950, "rRNA", "+"), ("ND6", 8950, 9400, "CDS", "+"),
                    ("ND3", 9395, 9600, "CDS", "+"), ("ND4L", 9600, 9850, "CDS", "+"),
                    ("ND1", 9850, 11200, "CDS", "+"), ("ND4", 11200, 12900, "CDS", "+"),
                    ("CYTB", 12900, 13900, "CDS", "+"), ("COX1", 13900, 16000, "CDS", "+")]]
    args = build_parser().parse_args(["plot", "--linear-layout", layout])
    resolve_linear_layout(args)
    fig = draw_publication_linear(args, reference, features, [], None, None)
    fig.canvas.draw()
    axis_text = {text.get_text() for text in fig.axes[0].texts}
    assert "1 bp" in axis_text and "0.001" not in axis_text
    labels = [text for text in fig.axes[0].texts if (text.get_gid() or "").startswith("annotation_label_")]
    assert sorted(t.get_text() for t in labels) == sorted(f["name"] for f in features)
    assert {t.get_fontsize() for t in labels} == {7}
    boxes = [t.get_window_extent() for t in labels]
    for i, box in enumerate(boxes):
        assert fig.bbox.contains(box.x0, box.y0) and fig.bbox.contains(box.x1, box.y1)
        assert not any(box.overlaps(other) for other in boxes[i + 1:])
    plt.close(fig)


def test_generic_repeat_does_not_acquire_terminal_cap_identity():
    args = build_parser().parse_args(["plot"])
    resolve_linear_layout(args)
    feature = dict(name="unrelated repeat", start=100, stop=500, type="repeat_region", strand="+", attributes={})
    fig = draw_publication_linear(args, Reference("mito", "ACGT" * 250), [feature], [], None, None)
    texts = {t.get_text() for t in fig.axes[0].texts}
    assert "Features" in texts
    assert "ITR" not in texts and "Cap" not in texts
    plt.close(fig)


@pytest.mark.parametrize("kind,name,attributes,present,absent", [
    ("repeat_region", "ITR core", {"rpt_type": "inverted"}, "ITR", "Cap"),
    ("misc_feature", "terminal cap", {}, "Cap", "ITR"),
])
def test_terminal_legend_only_names_supplied_feature_types(kind, name, attributes, present, absent):
    args = build_parser().parse_args(["plot"])
    resolve_linear_layout(args)
    feature = dict(name=name, start=0, stop=100, type=kind, strand="+", attributes=attributes)
    fig = draw_publication_linear(args, Reference("mito", "ACGT" * 1000), [feature], [], None, None)
    texts = {t.get_text() for t in fig.axes[0].texts + fig.texts}
    assert present in texts
    assert absent not in texts and "ITR / cap" not in texts
    plt.close(fig)


def test_multi_axes_svg_gradients_are_unique_and_read_width_floor_is_physical(tmp_path):
    segment = Segment("read", 0, 100, False, False, 60, [(0, 40), (2, 20), (0, 40)],
                      [(0, 40), (60, 100)], 0, 80, 80, 0, 0, False, "all", "")
    outline = read_outline(segment, 0, .7, 10, min_width=.25)
    assert np.min(np.abs(outline[:, 1])) == pytest.approx(.125)
    fig, axes = plt.subplots(2, figsize=(3.5, 2))
    for ax in axes:
        add_gradient_read(ax, segment, 0, .7, 10, ["#3a1d10", "#d3a878"])
        ax.set_xlim(0, 100)
        ax.set_ylim(-1, 1)
    args = build_parser().parse_args(["plot", "--fileform", "svg", "--no-timestamp", "-o", str(tmp_path / "multi")])
    save_linear_figure(fig, args, None, [])
    root = ET.parse(tmp_path / "multi.svg").getroot()
    ids = [node.get("id") for node in root.iter() if node.get("id")]
    assert len(ids) == len(set(ids))


def test_journal_minimum_does_not_erase_deletions():
    from matplotlib.path import Path
    from redwood.linear_redwood import GradientRead

    segment = Segment("read", 0, 100, False, False, 60, [(0, 40), (2, 20), (0, 40)],
                      [(0, 40), (60, 100)], 0, 80, 80, 0, 0, False, "all", "")
    args = build_parser().parse_args(["plot", "--publication-journal", "nature-communications"])
    resolve_linear_layout(args)
    fig = draw_publication_linear(args, Reference("mito", "A" * 100), [], [segment], None, None)
    read = fig.findobj(match=GradientRead)[0]
    # Normal aligned bases stay visibly thicker than the deletion, while
    # even the deletion's silhouette remains at least one physical point.
    coordinates = read.outline
    center = (coordinates[:, 1].min() + coordinates[:, 1].max()) / 2
    silhouette = Path(coordinates)
    assert silhouette.contains_point((20.5, center + .55))
    assert not silhouette.contains_point((50.5, center + .55))
    assert silhouette.contains_point((50.5, center + .499))
    assert silhouette.contains_point((50.5, center - .499))
    plt.close(fig)


@pytest.mark.parametrize("sign", [-1, 1])
def test_soft_clip_stack_has_disjoint_categories_and_cumulative_log_scale(sign):
    from matplotlib.colors import to_rgba
    from redwood.renderer import BARK_COLOR, REDWOOD_GRADIENT

    fig, ax = plt.subplots()
    short, long = np.array([3, 5, 0, 0]), np.array([6, 0, 8, 0])
    light, dark = add_clip_stack(ax, np.arange(4), np.ones(4), short, long, 15, sign, 6, 9)
    for i, (first, second) in enumerate(zip(light, dark)):
        assert first.get_y() == 15
        assert second.get_y() == pytest.approx(first.get_y() + first.get_height())
        assert first.get_height() == pytest.approx(sign * np.log1p(short[i]) / np.log(10) * 6)
        assert first.get_height() + second.get_height() == pytest.approx(sign * np.log1p(short[i] + long[i]) / np.log(10) * 6)
        a = sorted((first.get_y(), first.get_y() + first.get_height()))
        b = sorted((second.get_y(), second.get_y() + second.get_height()))
        assert max(a[0], b[0]) >= min(a[1], b[1]) - 1e-10
        assert first.get_facecolor() == to_rgba(REDWOOD_GRADIENT[-1])
        assert second.get_facecolor() == to_rgba(BARK_COLOR)
    plt.close(fig)


@pytest.mark.parametrize("layout", ["one-column", "two-column"])
def test_clip_length_legend_fits_without_overlapping_heading(layout):
    reference = Reference("mito", "ACGT" * 250)
    args = build_parser().parse_args(["plot", "--linear-layout", layout, "--linear-track", "clips"])
    resolve_linear_layout(args)
    profiles = {f"{kind}_{side}_clips": np.ones(reference.length, dtype=int)
                for kind in ("short", "long") for side in ("left", "right")}
    fig = draw_publication_linear(args, reference, [], [], {"profiles": profiles, "depth": {}, "summary": {}}, None)
    fig.canvas.draw()
    texts = [t for t in fig.axes[0].texts if t.get_text() in
             {"<100 bp", "≥100 bp", "Above: left / below: right", "max 50"}]
    assert len(texts) == 4
    boxes = [text.get_window_extent() for text in texts]
    for i, box in enumerate(boxes):
        assert fig.bbox.contains(box.x0, box.y0) and fig.bbox.contains(box.x1, box.y1)
        assert not any(box.overlaps(other) for other in boxes[i + 1:])
    plt.close(fig)
