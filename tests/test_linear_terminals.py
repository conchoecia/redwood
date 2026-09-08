from dataclasses import replace
from types import SimpleNamespace

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pytest

from redwood.linear_evidence import Reference, Segment, collect_evidence
from redwood.linear_terminals import _window_reads, draw_terminal_details


def segment(name, start, stop, left_clip=0, right_clip=0):
    cigar = ([(4, left_clip)] if left_clip else []) + [(0, stop - start)]
    cigar += [(4, right_clip)] if right_clip else []
    return Segment(name, start, stop, False, False, 60, cigar, [(start, stop)],
                   left_clip, left_clip + stop - start, left_clip + stop - start + right_clip,
                   left_clip, right_clip, False, "all", "")


def test_terminal_reads_are_unique_bounded_and_selected_before_ordering():
    reads = [segment(f"read{i}", i, 700 + 10 * i) for i in range(14)]
    reads += [segment("outside", 500, 1000), segment("touch_only", 700, 900)]
    reads += [replace(reads[4], supplementary=True)]
    chosen = _window_reads(reads, 0, 500)
    assert len(chosen) == 12
    assert [r[0] for r in chosen] == [f"read{i}" for i in range(2, 14)]
    assert len(next(r[1] for r in chosen if r[0] == "read4")) == 2


def test_terminal_figure_metadata_keeps_population_and_zoom_selection_separate():
    reference = Reference("mito", "ACGT" * 250)
    reads = [segment(f"read{i}", 0, 900 + i, 5, 0) for i in range(14)]
    evidence = collect_evidence(reads, reference, {"primary": 14, "supplementary": 0,
                                                 "secondary": 0, "hard_clipped": 0})
    before = evidence["profiles"]["starts"].copy()
    fig = draw_terminal_details(SimpleNamespace(width=3.5), reference, [], reads, evidence)
    try:
        metadata = fig._redwood_terminal_metadata
        assert metadata["layout"] == "stacked"
        assert [(p["start"], p["stop"]) for p in metadata["panels"]] == [(0, 500), (500, 1000)]
        assert all(len(p["read_ids"]) == 12 for p in metadata["panels"])
        assert all(p["profile_bin_size"] == 5 for p in metadata["panels"])
        assert any(t.get_text() == "Profiles: all 14 input reads" for ax in fig.axes for t in ax.texts)
        assert fig.axes[0].get_xlim() == (.5, 500.5)
        assert fig.axes[1].get_xlim() == (500.5, 1000.5)
        np.testing.assert_array_equal(before, evidence["profiles"]["starts"])
        assert metadata["panels"][0]["reads"][0]["segments"][0]["left_clip"] == 5
        fig.canvas.draw()
    finally:
        plt.close(fig)


def test_terminal_figure_short_reference_and_no_alignment_evidence():
    reference = Reference("short", "ACGT" * 5)
    fig = draw_terminal_details(SimpleNamespace(width=7.2), reference, [], [], None)
    try:
        metadata = fig._redwood_terminal_metadata
        assert metadata["layout"] == "side-by-side"
        assert all((p["first_base"], p["last_base"]) == (1, 20) for p in metadata["panels"])
        assert all(p["read_ids"] == [] for p in metadata["panels"])
        fig.canvas.draw()
    finally:
        plt.close(fig)


@pytest.mark.parametrize("journal,minimum", [("nature", .25), ("nature-communications", 1)])
def test_terminal_rules_and_indel_widths_follow_journal_minimum(journal, minimum):
    from redwood.linear_redwood import GradientRead

    reference = Reference("mito", "ACGT" * 250)
    read = segment("indel", 0, 1000, 5, 5)
    read.cigar = [(4, 5), (0, 20), (2, 40), (0, 940), (4, 5)]
    evidence = collect_evidence([read], reference, {"primary": 1, "supplementary": 0,
                                                  "secondary": 0, "hard_clipped": 0})
    fig = draw_terminal_details(SimpleNamespace(width=7.2, publication_journal=journal),
                                reference, [], [read], evidence)
    try:
        rules = fig.findobj(match=LineCollection)
        assert rules
        assert all(np.min(rule.get_linewidths()) >= minimum for rule in rules)
        for mesh in fig.findobj(match=GradientRead):
            y = mesh.outline[:, 1]
            center = (y.min() + y.max()) / 2
            assert min(abs(y - center)) * 2 >= minimum / 72 - 1e-12
        assert "all 1 input reads" in fig._redwood_caption
        assert fig._redwood_terminal_metadata["minimum_rule_pt"] == minimum
        fig.canvas.draw()
    finally:
        plt.close(fig)
