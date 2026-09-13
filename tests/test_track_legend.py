"""The ring key: a 90-degree cut-out of the track stack, one label per layer, drawn to the right of the map."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from redwood.renderer import (
    LEGEND_WIDTH,
    PLOT_LIMIT,
    add_spiral_read,
    add_track_legend,
    draw_circular_plot,
    track_legend_layers,
    units_per_point,
)
from redwood.multipass import RegularRead, pack_circular_reads


def test_layer_list_follows_what_is_drawn():
    keys = [k for k, _ in track_legend_layers(has_rnaseq=True, has_annotation=True, has_at=True, has_variants=True,
                                              has_numts=True, has_multipass=True, has_regular=True)]
    assert keys == ["regular", "multipass", "numts", "at", "variants", "genes", "trna", "rnaseq"]   # innermost first
    assert track_legend_layers(has_rnaseq=False, has_annotation=False, has_at=False, has_variants=False,
                               has_numts=False, has_multipass=False, has_regular=False) == []


def test_track_legend_labels_every_layer_and_stays_right_of_the_map():
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_aspect("equal"); ax.set_xlim(-PLOT_LIMIT, PLOT_LIMIT + LEGEND_WIDTH); ax.set_ylim(-PLOT_LIMIT, PLOT_LIMIT)
    layers = track_legend_layers(has_rnaseq=True, has_annotation=True, has_at=True, has_variants=True,
                                 has_numts=False, has_multipass=True, has_regular=True)
    n = add_track_legend(ax, layers)
    assert n == 7
    labels = [t.get_text() for t in ax.texts]
    for _, label in layers:
        assert label in labels
    # every label sits to the right of the map's own extent, and the artwork never reaches back into the map
    assert all(t.get_position()[0] > PLOT_LIMIT for t in ax.texts)
    for patch in ax.patches:
        xs = patch.get_path().vertices[:, 0]
        assert np.nanmin(xs) > PLOT_LIMIT
    plt.close(fig)


def test_circular_plot_widens_the_axes_only_with_the_legend():
    ref = "ACGT" * 500
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 5))
    draw_circular_plot(a1, length=2000, reference=ref, track_legend=True)
    draw_circular_plot(a2, length=2000, reference=ref, track_legend=False)
    assert a1.get_xlim() == (-PLOT_LIMIT, PLOT_LIMIT + LEGEND_WIDTH)
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
