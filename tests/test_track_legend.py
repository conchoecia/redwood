"""The ring key: a 90-degree cut-out of the track stack, one label per layer, in the bottom-right corner of the map."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from redwood.renderer import (
    KEY_CLEARANCE,
    PLOT_LIMIT,
    _key_center,
    add_spiral_read,
    add_track_legend,
    draw_circular_plot,
    track_legend_layers,
    track_legend_width,
    units_per_point,
)
from redwood.multipass import RegularRead, pack_circular_reads


def test_layer_list_follows_what_is_drawn():
    keys = [k for k, _ in track_legend_layers(has_rnaseq=True, has_annotation=True, has_at=True, has_variants=True,
                                              has_numts=True, has_multipass=True, has_regular=True)]
    assert keys == ["regular", "multipass", "numts", "at", "variants", "genes", "trna", "rnaseq"]   # innermost first
    assert track_legend_layers(has_rnaseq=False, has_annotation=False, has_at=False, has_variants=False,
                               has_numts=False, has_multipass=False, has_regular=False) == []


def _key_axes(figsize=(6.3, 5.8)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect("equal"); ax.set_ylim(-PLOT_LIMIT, PLOT_LIMIT)
    ax.set_xlim(-PLOT_LIMIT, PLOT_LIMIT + track_legend_width(ax))
    return fig, ax


def _data_points(ax, artist):
    """Rendered outline of an artist in data coordinates (curves flattened, so Bezier control points do not count)."""
    if hasattr(artist, "get_xydata"):
        return artist.get_xydata()
    to_data = artist.get_transform() - ax.transData
    return np.concatenate(artist.get_path().to_polygons(transform=to_data, closed_only=False))


def test_track_legend_sits_in_the_corner_clear_of_the_map():
    fig, ax = _key_axes()
    layers = track_legend_layers(has_rnaseq=True, has_annotation=True, has_at=True, has_variants=True,
                                 has_numts=False, has_multipass=True, has_regular=True)
    assert add_track_legend(ax, layers) == 7
    labels = [t.get_text() for t in ax.texts]
    assert all(label in labels for _, label in layers) and "ring key" in labels
    assert not any("cent" in t.lower() for t in labels)                  # no caption for the key's center
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    for artist in list(ax.patches) + list(ax.lines):
        xy = _data_points(ax, artist)
        assert np.hypot(xy[:, 0], xy[:, 1]).min() >= KEY_CLEARANCE - 1e-3  # never reaches into the rings or their labels
        assert x0 <= xy[:, 0].min() and xy[:, 0].max() <= x1 and y0 <= xy[:, 1].min() and xy[:, 1].max() <= y1
        assert xy[:, 0].min() > 0 and xy[:, 1].max() < 0                  # bottom-right quadrant
    for t in ax.texts:
        assert 3.6 <= t.get_fontsize() <= 5.0 and t.get_position()[0] > 0 and t.get_position()[1] < 0
    plt.close(fig)


def test_multipass_key_sample_spirals_inward():
    fig, ax = _key_axes()
    add_track_legend(ax, [("multipass", "multi-pass reads")])
    cx, cy = _key_center()
    assert len(ax.lines) == 3
    for line in ax.lines:
        xy = line.get_xydata()
        ang = np.degrees(np.arctan2(xy[:, 1] - cy, xy[:, 0] - cx))
        r = np.hypot(xy[:, 0] - cx, xy[:, 1] - cy)
        assert np.all(np.diff(ang) <= 1e-9)                   # travels clockwise, like the map
        assert np.all(np.diff(r) <= 1e-9)                     # and only ever steps inward
    plt.close(fig)


def test_circular_plot_widens_the_axes_only_with_the_legend():
    ref = "ACGT" * 500
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 5))
    draw_circular_plot(a1, length=2000, reference=ref, track_legend=True)
    draw_circular_plot(a2, length=2000, reference=ref, track_legend=False)
    assert PLOT_LIMIT < a1.get_xlim()[1] < PLOT_LIMIT + 0.6          # a narrow label column, not a second panel
    assert a2.get_xlim() == (-PLOT_LIMIT, PLOT_LIMIT)
    assert any(t.get_text() == "AT content" for t in a1.texts) and not any(t.get_text() == "AT content" for t in a2.texts)
    plt.close(fig)


def test_spiral_linewidth_never_exceeds_the_rung_spacing():
    fig, ax = plt.subplots(figsize=(2.5, 2.5))          # a small panel, where a fixed 0.9 pt line would merge the rungs
    ax.set_aspect("equal"); ax.set_xlim(-PLOT_LIMIT, PLOT_LIMIT); ax.set_ylim(-PLOT_LIMIT, PLOT_LIMIT)
    rung = 0.0087
    add_spiral_read(ax, 0, 3.0, 1000, 0.8, rung, 0.12, "#d1322b")
    lw = ax.lines[-1].get_linewidth()
    assert lw <= 0.6 * rung / units_per_point(ax) + 1e-9 and lw >= 0.25
    plt.close(fig)


def test_full_length_reads_pack_in_start_order():
    reads = [RegularRead(start, 1000) for start in (700, 100, 400)]
    placed, n_rungs = pack_circular_reads(reads, 1000)
    assert n_rungs == 3 and [rd.start for rd, _ in placed] == [100, 400, 700]
