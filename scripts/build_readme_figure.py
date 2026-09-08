#!/usr/bin/env python
"""Build the light/dark README preview figure from the package renderer."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/redwood-matplotlib")

import matplotlib.pyplot as plt

from redwood.renderer import draw_dataset_panel
from redwood.cli import build_parser
from redwood.linear_evidence import (
    collect_evidence, load_annotations, load_reference, read_segments, select_display_segments,
)
from redwood.linear_publication import draw_publication_linear, publication_font, resolve_linear_layout


DATASETS = [
    ("human", "Human", "Homo sapiens", 80),
    ("mouse", "Mouse", "Mus musculus", 80),
    ("drosophila", "Fly", "Drosophila melanogaster", 80),
    ("sponge", "Sponge", "Ephydatia muelleri", 80),
]


def write_grid(datasets_dir: Path, output: Path, dark: bool, rnaseq_style: str,
               overview: bool = False) -> None:
    bg = "#0d1117" if dark else "#ffffff"
    edge = "#303946" if dark else "#d8dee8"
    datasets = [DATASETS[0], DATASETS[3]] if overview else DATASETS
    fig, axes = plt.subplots(1 if overview else 2, 2,
                             figsize=(10, 5.4 if overview else 11), dpi=170)
    fig.patch.set_facecolor(bg)
    for ax, (dataset, label, species, max_reads) in zip(axes.flat, datasets):
        ax.set_facecolor(bg)
        for spine in ax.spines.values():
            spine.set_color(edge)
            spine.set_linewidth(1.0)
        draw_dataset_panel(ax, datasets_dir / dataset, label, species, max_reads, dark, rnaseq_style)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(left=0.035, right=0.965, bottom=0.055, top=0.95, wspace=0.10, hspace=0.30)
    fig.savefig(output, facecolor=bg)
    plt.close(fig)


def write_linear(datasets_dir: Path, output: Path, dark: bool) -> None:
    """Use the publication renderer directly, with a small README title margin."""
    dataset = datasets_dir / "hydra"
    manifest = json.loads((dataset / "manifest.json").read_text())
    args = build_parser().parse_args([
        "plot", "--topology", "linear", "--linear-layout", "two-column",
        "--max-reads", "30", "--terminal-window", "100",
        "--linear-track", "depth", "--linear-track", "ends",
        "--linear-track", "clips", "--no-timestamp", "-T",
    ] + (["--dark"] if dark else []))
    resolve_linear_layout(args)
    reference = load_reference(dataset / manifest["reference"])
    features = load_annotations(dataset / manifest["annotation"], reference)
    segments, counts = read_segments(dataset / manifest["bam"], reference, None)
    evidence = collect_evidence(segments, reference, counts, args.terminal_window,
                                args.junction_window, args.clip_threshold)
    selected = select_display_segments(segments, reference.length, args.query, args.sort,
                                       args.max_reads, args.linear_read_selection, args.terminal_window)
    with plt.rc_context({"font.family": publication_font()}):
        fig = draw_publication_linear(args, reference, features, selected, evidence, None)
        width, height = fig.get_size_inches()
        margin = .44
        fig.set_size_inches(width, height + margin)
        # Keep the plot's physical font and line sizes while making room for
        # the README heading above it. Figure-relative gutter labels move too.
        factor = height / (height + margin)
        for ax in fig.axes:
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 * factor, box.width, box.height * factor])
        for label in fig.texts:
            x, y = label.get_position()
            label.set_position((x, y * factor))
        fg = "#eef4fb" if dark else "#20252b"
        fig.text(.5, 1 - .11 / (height + margin), "Hydra", ha="center", va="top",
                 fontsize=12, weight="bold", color=fg)
        fig.text(.5, 1 - .29 / (height + margin), "Hydra oligactis · linear mitochondrial genome",
                 ha="center", va="top", fontsize=7.5, style="italic", color=fg)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=240, facecolor=fig.get_facecolor())
        plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets-dir", type=Path, default=Path("examples/datasets"))
    parser.add_argument("--outdir", type=Path, default=Path("docs/assets"))
    parser.add_argument("--only", choices=["all", "overview", "linear", "gallery"], default="all",
                        help="Build all previews (default), or one part while adjusting the README.")
    parser.add_argument(
        "--rnaseq-style",
        choices=["coverage", "strand"],
        default="coverage",
        help="Render RNA-seq as a plain coverage histogram or color it by alignment strand.",
    )
    args = parser.parse_args(argv)

    for dark in (False, True):
        theme = "dark" if dark else "light"
        if args.only in {"all", "gallery"}:
            write_grid(args.datasets_dir, args.outdir / f"redwood-grid-{theme}.png", dark, args.rnaseq_style)
        if args.only in {"all", "overview"}:
            write_grid(args.datasets_dir, args.outdir / f"redwood-circular-{theme}.png", dark, args.rnaseq_style, overview=True)
        if args.only in {"all", "linear"}:
            write_linear(args.datasets_dir, args.outdir / f"redwood-linear-{theme}.png", dark)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
