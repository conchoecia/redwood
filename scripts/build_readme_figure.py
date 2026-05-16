#!/usr/bin/env python
"""Build the light/dark README preview figure from the package renderer."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/redwood-matplotlib")

import matplotlib.pyplot as plt

from redwood.renderer import draw_dataset_panel


DATASETS = [
    ("human", "Human", "Homo sapiens", 24),
    ("mouse", "Mouse", "Mus musculus", 36),
    ("drosophila", "Fly", "Drosophila melanogaster", 44),
    ("sponge", "Sponge", "Ephydatia muelleri", 56),
]


def write_grid(datasets_dir: Path, output: Path, dark: bool, rnaseq_style: str) -> None:
    bg = "#0d1117" if dark else "#ffffff"
    edge = "#303946" if dark else "#d8dee8"
    fig, axes = plt.subplots(2, 2, figsize=(10, 11), dpi=170)
    fig.patch.set_facecolor(bg)
    for ax, (dataset, label, species, max_reads) in zip(axes.flat, DATASETS):
        ax.set_facecolor(bg)
        for spine in ax.spines.values():
            spine.set_color(edge)
            spine.set_linewidth(1.0)
        draw_dataset_panel(ax, datasets_dir / dataset, label, species, max_reads, dark, rnaseq_style)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(left=0.035, right=0.965, bottom=0.055, top=0.95, wspace=0.10, hspace=0.30)
    fig.savefig(output, facecolor=bg)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets-dir", type=Path, default=Path("examples/datasets"))
    parser.add_argument("--outdir", type=Path, default=Path("docs/assets"))
    parser.add_argument(
        "--rnaseq-style",
        choices=["coverage", "strand"],
        default="coverage",
        help="Render RNA-seq as a plain coverage histogram or color it by alignment strand.",
    )
    args = parser.parse_args(argv)

    write_grid(args.datasets_dir, args.outdir / "redwood-grid-light.png", dark=False, rnaseq_style=args.rnaseq_style)
    write_grid(args.datasets_dir, args.outdir / "redwood-grid-dark.png", dark=True, rnaseq_style=args.rnaseq_style)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
