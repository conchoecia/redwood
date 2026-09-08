#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import contextlib
import io
import os
from pathlib import Path

from .renderer import run_plot
from .workflow import (
    map_long,
    map_rnaseq,
    prepare_reference,
    run_end_to_end,
    write_metrics,
)


class FullPaths(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def build_parser():
    parser = argparse.ArgumentParser(
        prog="redwood",
        description="Plot circular or linear genome read, annotation, and depth tracks.",
    )
    subparsers = parser.add_subparsers(dest="command")

    parser_plot = subparsers.add_parser(
        "plot",
        help="make a redwood genome plot",
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
    add_topology_argument(parser_plot)
    add_linear_arguments(parser_plot)
    parser_plot.add_argument("--reprocess-ont", action="store_true",
                             help="Treat the main BAM as ONT: extract mapped molecules, trim adapters, and remap before plotting.")
    parser_plot.add_argument("--ont-workdir", type=Path,
                             help="New directory for ONT reads, trimming reports and remapped BAM (default: OUTPUT.ont-preprocess).")
    add_ont_arguments(parser_plot)
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
    parser_plot.add_argument(
        "--mito-fasta",
        "--reference-fasta",
        dest="mito_fasta",
        metavar="fasta",
        action=FullPaths,
        help="Mitochondrial FASTA used for sequence-derived tracks.",
    )
    parser_plot.add_argument("-I", "--interlace", action="store_true", default=False)
    parser_plot.add_argument("-i", "--invert", action="store_true", default=False)
    parser_plot.add_argument("-L", "--log", action="store_true", default=False)
    parser_plot.add_argument("-M", "--main-bam", dest="main_bam", metavar="mainbam", action=FullPaths)
    parser_plot.add_argument("--max-reads", type=int,
                             help="Displayed reads (default: 30 for one-column linear, otherwise 80).")
    parser_plot.add_argument(
        "--max-internal-gap",
        dest="max_internal_gap",
        type=int,
        default=50,
        help="Largest internal alignment gap (bp) a multi-pass read may "
        "contain; a larger gap is treated as adapter/chimeric junk.",
    )
    parser_plot.add_argument(
        "--min-pass-fraction",
        dest="min_pass_fraction",
        type=float,
        default=1.0,
        help="Minimum alignment span, in circle lengths, for a read to be "
        "drawn as a multi-pass spiral; 1.0 = any read covering the full circle.",
    )
    parser_plot.add_argument(
        "--wrap-ramp",
        dest="wrap_ramp",
        type=float,
        default=0.12,
        help="Fraction of each spiral turn used for the diagonal step-down "
        "into the next rung.",
    )
    parser_plot.add_argument(
        "--no-multipass",
        dest="no_multipass",
        action="store_true",
        help="Disable multi-pass spiral detection; draw every read single-pass.",
    )
    parser_plot.add_argument(
        "--min-indel",
        dest="min_indel",
        type=int,
        default=10,
        help="Smallest insertion/deletion (bp) drawn as an indel on a read "
        "arc; shorter ones render as match. Use 1 to draw every indel.",
    )
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
        default=None,
        nargs="+",
        help="Pandas query clauses for displayed linear read rows (e.g. 'ALNLEN >= 10000' "
        "'MAPLEN <= reflength'). Depth and evidence always use all input reads.",
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
        help="Circular read order; linear length ranking (POS uses ALNLEN). Linear rows always follow reference position.",
    )
    parser_plot.add_argument("--ticks", type=int, nargs="+", default=[0, 10, 100, 1000])
    parser_plot.add_argument("--title")
    parser_plot.add_argument("--subtitle")
    parser_plot.add_argument("--dark", action="store_true", help="Render using the dark plot theme.")
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
    parser_plot.set_defaults(func=run_plot)

    parser_advanced = subparsers.add_parser(
        "advanced",
        help="lower-level workflow steps for debugging and custom pipelines",
    )
    advanced_subparsers = parser_advanced.add_subparsers(dest="advanced_command")
    advanced_subparsers.required = True

    parser_prepare = advanced_subparsers.add_parser(
        "prepare-reference",
        help="write derived mitochondrial and RNA-seq bait references",
    )
    parser_prepare.add_argument("--mito-fasta", required=True, type=Path)
    parser_prepare.add_argument("--nuclear-fasta", type=Path)
    parser_prepare.add_argument("--outdir", required=True, type=Path)
    add_topology_argument(parser_prepare)
    parser_prepare.add_argument(
        "--exclude-token",
        action="append",
        default=[],
        help="Additional case-insensitive nuclear FASTA header token to exclude from RNA-seq bait references.",
    )
    parser_prepare.set_defaults(func=prepare_reference)

    parser_long = advanced_subparsers.add_parser(
        "map-long",
        help="map long reads using the requested topology and select circular plot reads",
    )
    parser_long.add_argument("--mito-fasta", required=True, type=Path)
    parser_long.add_argument("--long-reads", required=True, type=Path, nargs="+")
    parser_long.add_argument("--outdir", required=True, type=Path)
    add_topology_argument(parser_long)
    parser_long.add_argument("--output-bam", type=Path)
    parser_long.add_argument(
        "--preset",
        choices=["map-ont", "map-pb", "map-hifi", "asm5", "asm10", "asm20"],
        default="map-ont",
        help="minimap2 preset for long-read mapping.",
    )
    parser_long.add_argument(
        "--copies",
        default="auto",
        help="Tandem copies in the mapping reference: an integer, or 'auto' "
        "to size from the longest read. Reads that circle the genome "
        "multiple times need a reference long enough to hold them in one "
        "continuous alignment.",
    )
    parser_long.add_argument("--target-depth", type=float, default=100.0)
    parser_long.add_argument("--min-span-fraction", type=float, default=0.25)
    parser_long.add_argument("--dry-run", action="store_true")
    add_ont_arguments(parser_long, mapping=True)
    parser_long.set_defaults(func=map_long)

    from .ont import reprocess_existing

    parser_ont = advanced_subparsers.add_parser(
        "reprocess-ont", help="trim and remap ONT molecules already mapped to the target genome")
    parser_ont.add_argument("--mito-fasta", required=True, type=Path)
    parser_ont.add_argument("--main-bam", required=True, type=Path)
    parser_ont.add_argument("--outdir", required=True, type=Path)
    add_topology_argument(parser_ont)
    add_ont_arguments(parser_ont)
    parser_ont.set_defaults(func=reprocess_existing)

    parser_rna = advanced_subparsers.add_parser(
        "map-rnaseq",
        help="map RNA-seq to nuclear bait plus mitochondrion and keep mitochondrial alignments",
    )
    parser_rna.add_argument("--mito-fasta", required=True, type=Path)
    parser_rna.add_argument("--nuclear-fasta", required=True, type=Path)
    parser_rna.add_argument("--rnaseq-reads", required=True, type=Path, nargs="+")
    parser_rna.add_argument("--outdir", required=True, type=Path)
    parser_rna.add_argument("--output-bam", type=Path)
    parser_rna.add_argument("--preset", default="sr", help="minimap2 preset for RNA-seq mapping.")
    parser_rna.add_argument("--exclude-token", action="append", default=[])
    parser_rna.add_argument("--dry-run", action="store_true")
    parser_rna.set_defaults(func=map_rnaseq)

    parser_metrics = subparsers.add_parser(
        "metrics",
        help="summarize long-read and RNA-seq mitochondrial support metrics",
    )
    parser_metrics.add_argument("--mito-fasta", required=True, type=Path)
    parser_metrics.add_argument("--long-bam", type=Path)
    parser_metrics.add_argument("--rnaseq-bam", type=Path)
    parser_metrics.add_argument("--output", required=True, type=Path)
    add_topology_argument(parser_metrics)
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
    add_topology_argument(parser_run)
    add_linear_arguments(parser_run)
    add_ont_arguments(parser_run, mapping=True)
    parser_run.add_argument("--long-read-preset", default="map-ont")
    parser_run.add_argument("--rnaseq-preset", default="sr")
    parser_run.add_argument("--long-read-depth", type=float, default=100.0)
    parser_run.add_argument("--max-reads", type=int,
                            help="Displayed reads (default: 30 for one-column linear, otherwise 80).")
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


def add_ont_arguments(parser, mapping=False):
    if mapping:
        parser.add_argument("--no-ont-trim", action="store_true",
                            help="Skip automatic adapter trimming/remapping with the map-ont preset.")
    parser.add_argument("--ont-adapter-5p", action="append",
                        help="5-prime adapter in original read orientation; repeat to replace the default ONT ligation motifs.")
    parser.add_argument("--ont-adapter-3p", action="append",
                        help="3-prime adapter in original read orientation; repeat to replace the default ONT ligation motifs.")
    parser.add_argument("--ont-overlap", type=int, default=12,
                        help="Minimum adapter overlap (default: 12 bp).")
    parser.add_argument("--ont-error-rate", type=float, default=0.1,
                        help="Cutadapt maximum substitution/indel error rate (default: 0.1).")
    parser.add_argument("--ont-end-window", type=int, default=100,
                        help="Search only this many bases at each read end (default: 100; capped at half the read length).")
    parser.add_argument("--ont-threads", type=int, default=4,
                        help="Threads for ONT trimming and remapping (default: 4).")
    parser.add_argument("--ont-source-reads", type=Path, nargs="+",
                        help="Original ONT FASTQ(s), required if complete read sequence/qualities cannot be recovered from BAM.")


def add_topology_argument(parser):
    parser.add_argument("--topology", choices=["circular", "linear"], default="circular",
                        help="Genome topology. Linear mode never doubles or wraps the reference.")


def add_linear_arguments(parser):
    parser.add_argument("--linear-layout", choices=["one-column", "two-column", "legacy"], default="two-column",
                        help="Redwood linear layout at publication size (default: two-column); legacy keeps the 13-inch layout.")
    parser.add_argument("--publication-journal", choices=["nature", "nature-communications"], default="nature",
                        help="Column widths and minimum line weights for publication layouts (default: nature).")
    parser.add_argument("--panel-label", help="Optional publication panel letter.")
    parser.add_argument("--terminal-details", action="store_true",
                        help="Also export a companion figure with expanded terminal coordinates.")
    parser.add_argument("--variant-sites", type=Path,
                        help="Also export an allele panel for a TSV of 1-based sites (position or pos column).")
    parser.add_argument("--variant-min-base-quality", type=int, default=20,
                        help="Minimum base quality for the optional allele panel (default: 20).")
    parser.add_argument("--linear-read-selection", choices=["terminal-balanced", "longest"],
                        default="terminal-balanced",
                        help="Linear read sample: longest reads balanced across left/right/both termini (default), "
                             "or longest overall. Selection precedes coordinate sorting.")
    parser.add_argument("--linear-style", choices=["redwood", "diagnostic"], default="redwood",
                        help="Linear figure style: Redwood's existing tracks unrolled (default), or diagnostic panels.")
    parser.add_argument("--linear-track", choices=["depth", "ends", "clips", "none"], action="append",
                        help="Production evidence bands (default: depth). Repeat for depth/ends/clips, or use none.")
    parser.add_argument("--read-color", choices=["wood", "class"], default="wood",
                        help="Production linear reads: original wood gradient (default) or supplied class colors.")
    parser.add_argument("--show-terminal-sequences", action="store_true",
                        help="Include the first/last 30 bases beneath the production linear figure.")
    parser.add_argument("--rnaseq-label", default="RNA depth",
                        help="RNA track label in production linear figures (default: RNA depth).")
    parser.add_argument("--read-classes", type=Path,
                        help="Linear plots: TSV with read and class columns, optional locus/nuclear_locus.")
    parser.add_argument("--terminal-window", type=int, default=30,
                        help="Linear read selection and evidence: terminal window in bp (default: 30).")
    parser.add_argument("--junction-window", type=int, default=300,
                        help="Linear evidence: distance from a terminus for split joins (default: 300).")
    parser.add_argument("--clip-threshold", type=int, default=100,
                        help="Linear evidence: minimum long soft-clip length in bp (default: 100).")
    parser.add_argument("--bin-size", type=int, default=25,
                        help="Linear plots: start/end and soft-clip histogram bin size in bp (default: 25).")
    parser.add_argument("--depth-scale", choices=["linear", "log"], default="linear",
                        help="Linear plots: depth axis scale (log uses log(1 + depth), retaining zeros).")
    parser.add_argument("--hide-evidence", action="store_true",
                        help="Linear plots: omit start/end and clip panels; still export evidence tables.")
    parser.add_argument("--width", type=float,
                        help="Override linear width in inches; otherwise use the selected publication layout.")


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 2
    if args.command != "plot" or getattr(args, "verbose", False) or getattr(args, "topology", "circular") == "linear":
        args.func(args)
    else:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
