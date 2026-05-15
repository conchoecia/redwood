#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import contextlib
import io
import os

from .plot import run


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
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 2
    if args.verbose:
        args.func(args)
    else:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
