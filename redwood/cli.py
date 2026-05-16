#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import contextlib
import io
import os
from pathlib import Path

from .plot import run
from .workflow import (
    run_end_to_end,
    write_metrics,
)


class FullPaths(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def build_parser():
    parser = argparse.ArgumentParser(
        prog="redwood",
        description="Plot circular genome read, annotation, and depth tracks.",
    )
    subparsers = parser.add_subparsers(dest="command")

    parser_plot = subparsers.add_parser(
        "plot",
        help="make a redwood circular genome plot",
    )
    parser_plot.add_argument(
        "-d",
        "--doubled",
        dest="doubled",
        choices=["main", "rnaseq"],
        default=[],
        nargs="+",
        help="Input BAMs mapped to a doubled circular reference. Accepts main, rnaseq, or both.",
    )
    parser_plot.add_argument("--dpi", metavar="dpi", default=600, type=int)
    parser_plot.add_argument(
        "--fileform",
        dest="fileform",
        metavar="STRING",
        choices=[
            "png",
            "pdf",
            "eps",
            "jpeg",
            "jpg",
            "pgf",
            "ps",
            "raw",
            "rgba",
            "svg",
            "svgz",
            "tif",
            "tiff",
        ],
        default=["png"],
        nargs="+",
    )
    parser_plot.add_argument("--gff", metavar="gff", action=FullPaths)
    parser_plot.add_argument("-I", "--interlace", action="store_true", default=False)
    parser_plot.add_argument("-i", "--invert", action="store_true", default=False)
    parser_plot.add_argument("-L", "--log", action="store_true", default=False)
    parser_plot.add_argument("-M", "--main-bam", dest="main_bam", metavar="mainbam", action=FullPaths)
    parser_plot.add_argument("--no-timestamp", dest="no_timestamp", action="store_true")
    parser_plot.add_argument(
        "-o",
        "--output-base-name",
        dest="BASENAME",
        help="Base name for output files. Defaults to redwood.",
    )
    parser_plot.add_argument(
        "--query",
        dest="query",
        default=["ALNLEN >= 10000", "MAPLEN < reflength"],
        nargs="+",
        help="Pandas query clauses used to filter long-read BAM rows.",
    )
    parser_plot.add_argument("-R", "--rnaseq-bam", dest="rnaseq_bam", metavar="rnabam", action=FullPaths)
    parser_plot.add_argument(
        "--small-start",
        dest="small_start",
        choices=["inside", "outside"],
        default="inside",
    )
    parser_plot.add_argument(
        "--sort",
        dest="sort",
        choices=["ALNLEN", "TRULEN", "MAPLEN", "POS"],
        default="ALNLEN",
    )
    parser_plot.add_argument("--ticks", type=int, nargs="+", default=[0, 10, 100, 1000])
    parser_plot.add_argument(
        "--extra-track",
        dest="extra_tracks",
        choices=["at", "gc", "rnaseq-strand", "metrics"],
        default=[],
        nargs="+",
        help=(
            "Declare optional tracks for newer plot styles. The legacy plotter "
            "currently uses the default read, annotation, and RNA-seq depth tracks."
        ),
    )
    parser_plot.add_argument(
        "-T",
        "--transparent",
        action="store_false",
        default=True,
        help="Use an opaque background. Default output background is transparent.",
    )
    parser_plot.add_argument(
        "--verbose",
        action="store_true",
        help="Show progress/debug output from the plotter.",
    )
    parser_plot.set_defaults(func=run)

    parser_metrics = subparsers.add_parser(
        "metrics",
        help="summarize long-read and RNA-seq mitochondrial support metrics",
    )
    parser_metrics.add_argument("--mito-fasta", required=True, type=Path)
    parser_metrics.add_argument("--long-bam", type=Path)
    parser_metrics.add_argument("--rnaseq-bam", type=Path)
    parser_metrics.add_argument("--output", required=True, type=Path)
    parser_metrics.set_defaults(func=write_metrics)

    parser_run = subparsers.add_parser(
        "run",
        help="run an end-to-end local redwood workflow from references and reads",
    )
    parser_run.add_argument("--mito-fasta", required=True, type=Path)
    parser_run.add_argument(
        "--nuclear-fasta",
        type=Path,
        help="Whole-genome FASTA used as RNA-seq bait; mitochondrial-looking contigs are excluded.",
    )
    parser_run.add_argument("--gff", type=Path, help="Optional GFF3 annotation for the mitochondrial genome.")
    parser_run.add_argument("--long-reads", type=Path, nargs="+")
    parser_run.add_argument("--rnaseq-reads", type=Path, nargs="+")
    parser_run.add_argument("--outdir", required=True, type=Path)
    parser_run.add_argument("--long-read-preset", default="map-ont")
    parser_run.add_argument("--rnaseq-preset", default="sr")
    parser_run.add_argument("--long-read-depth", type=float, default=50.0)
    parser_run.add_argument("--min-span-fraction", type=float, default=0.25)
    parser_run.add_argument("--exclude-token", action="append", default=[])
    parser_run.add_argument("--plot-name", default="redwood")
    parser_run.add_argument("--skip-plot", action="store_true")
    parser_run.add_argument("--dry-run", action="store_true")
    parser_run.add_argument("--dpi", default=600, type=int)
    parser_run.add_argument("--fileform", default=["png"], nargs="+")
    parser_run.add_argument("--ticks", type=int, nargs="+", default=[0, 10, 100, 1000])
    parser_run.add_argument(
        "-T",
        "--transparent",
        action="store_false",
        default=True,
        help="Use an opaque background. Default output background is transparent.",
    )
    parser_run.set_defaults(func=run_end_to_end)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 2
    if args.command != "plot" or getattr(args, "verbose", False):
        args.func(args)
    else:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
