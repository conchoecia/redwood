#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import contextlib
import io
import os
from pathlib import Path

from .renderer import run_plot
from .numts import run_numts
from .composite import STYLES, run_composite
from .variants import run_variants
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
        "--trna-labels",
        dest="trna_labels",
        choices=["letter", "name", "none"],
        default="letter",
        help="How to label tRNAs: one-letter amino-acid code (default), the full GFF name, or no label.",
    )
    parser_plot.add_argument(
        "--no-feature-labels",
        dest="no_feature_labels",
        action="store_true",
        help="Do not label CDS/rRNA/gene features.",
    )
    parser_plot.add_argument(
        "--read-mismatches", dest="read_mismatches", choices=["shared", "all", "none"], default="shared",
        help="IGV-style marks on the read rings (needs --mito-fasta): 'shared' marks mismatches only at columns "
             "where the read population disagrees with the reference above --min-minor-frac, 'all' marks every "
             "mismatch, 'none' disables them. Insertions/deletions >= --min-indel are always marked.",
    )
    parser_plot.add_argument("--numt-loci", dest="numt_loci", action=FullPaths,
                             help="redwood numts loci TSV: draws a ring of the mitogenome intervals present as NUMTs, colored by identity.")
    parser_plot.add_argument("--circular-read-classes", dest="circular_read_classes", action=FullPaths,
                             help="redwood numts read-class TSV: color read arcs by class (NUMT junction reads red, chimeras purple, ...).")
    parser_plot.add_argument("--no-variant-ring", dest="no_variant_ring", action="store_true",
                             help="Do not draw the per-column disagreement ring inside the annotation.")
    parser_plot.add_argument("--no-track-legend", dest="no_track_legend", action="store_true",
                             help="Omit the ring key (90-degree cut-out of the track stack drawn to the right of the map).")
    parser_plot.add_argument("--variant-table", dest="variant_table", action=FullPaths,
                             help="Use this redwood variants TSV (e.g. from the whole read set) instead of computing "
                                  "column disagreement from --main-bam.")
    parser_plot.add_argument("--min-minor-frac", dest="min_minor_frac", type=float, default=0.05,
                             help="Minor-allele fraction that flags a column (default 0.05; scaled up to 3x the "
                                  "table-wide median for noisy reads, see `redwood variants --noise-multiplier`).")
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

    parser_variants = subparsers.add_parser(
        "variants",
        help="per-column table of read mismatches, insertions and deletions against the mitogenome",
    )
    parser_variants.add_argument("--mito-fasta", required=True, type=Path)
    parser_variants.add_argument("--bam", required=True, type=Path,
                                 help="Reads mapped to the mitogenome (single-copy or the doubled reference).")
    parser_variants.add_argument("--output", required=True, type=Path, help="TSV of flagged columns.")
    parser_variants.add_argument("--summary", type=Path, help="Optional JSON summary.")
    parser_variants.add_argument("--contig", help="Restrict to one BAM contig (default: all).")
    parser_variants.add_argument("--all-columns", action="store_true", help="Write every column, not only flagged ones.")
    parser_variants.add_argument("--min-base-quality", type=int, default=20)
    parser_variants.add_argument("--min-depth", type=int, default=5)
    parser_variants.add_argument("--min-minor-frac", type=float, default=0.05,
                                 help="Flag a column as 'minor' / 'insertion' when the fraction reaches this (default 0.05).")
    parser_variants.add_argument("--noise-multiplier", type=float, default=3.0,
                                 help="Flag thresholds = max(--min-minor-frac, this x the table-wide median fraction), so "
                                      "noisy CLR/ONT reads flag only columns above their error level (0 = fixed threshold).")
    parser_variants.set_defaults(func=run_variants)

    parser_numts = subparsers.add_parser(
        "numts",
        help="catalog NUMTs in a nuclear assembly, classify reads as mitochondrial or NUMT-derived, draw NUMT figures",
    )
    parser_numts.add_argument("--mito-fasta", required=True, type=Path)
    parser_numts.add_argument("--nuclear-fasta", required=True, type=Path, help="Nuclear assembly (exclude any mitochondrial scaffold first).")
    parser_numts.add_argument("--outdir", required=True, type=Path)
    parser_numts.add_argument("--gff", type=Path, help="Mitogenome annotation (drawn under the catalog panel).")
    parser_numts.add_argument("--long-reads", type=Path, nargs="+", help="Long reads to map to nuclear + mito and classify.")
    parser_numts.add_argument("--long-read-preset", default="map-hifi")
    parser_numts.add_argument("--bam", type=Path, help="Existing BAM of reads vs nuclear + mito (minimap2 -Y --secondary=no) instead of --long-reads.")
    parser_numts.add_argument("--mito-bam", type=Path, help="Reads mapped to the mitogenome alone; a copy tagged PO:Z=<class> is written.")
    parser_numts.add_argument("--blastn", action="store_true", help="Also run blastn -task dc-megablast for short/diverged fragments (needs BLAST+).")
    parser_numts.add_argument("--merge", type=int, default=3000, help="Merge hits within this distance into one locus (default 3000).")
    parser_numts.add_argument("--pad", type=int, default=2000, help="A nuclear read segment within this distance of a locus counts as 'at the locus'.")
    parser_numts.add_argument("--min-flank", dest="min_flank", type=int, default=500,
                              help="A mito+nuclear read counts as a NUMT junction read only with this much nuclear sequence outside the locus.")
    parser_numts.add_argument("--flank", type=int, default=1000, help="Junction support: alignment must extend this far on both sides of a locus end.")
    parser_numts.add_argument("--threads", type=int, default=8)
    parser_numts.add_argument("--dpi", type=int, default=200)
    parser_numts.add_argument("--figure", action="store_true", help="Also draw the composite figure: circular map, landscape, catalog.")
    parser_numts.add_argument("--figure-name", default="mitogenome_numts")
    parser_numts.add_argument("--figure-bam", type=Path, help="Long-read BAM for the circular panel (e.g. redwood long_reads.redwood.bam).")
    parser_numts.add_argument("--figure-rnaseq-bam", type=Path)
    parser_numts.add_argument("--figure-variant-table", type=Path)
    parser_numts.set_defaults(func=run_numts)

    parser_comp = subparsers.add_parser(
        "composite",
        help="journal-style figure: circular map + NUMT panels, as a figure or a full page with the legend",
        description="Draw the mitogenome map with the NUMT panels laid out to a journal style. Styles set the figure width, text "
                    "and panel-label sizes, fonts, page geometry and legend format; every setting can be overridden.",
    )
    gi = parser_comp.add_argument_group("inputs")
    gi.add_argument("--mito-fasta", required=True, type=Path)
    gi.add_argument("--numt-loci", required=True, type=Path, help="numts.loci.tsv from redwood numts.")
    gi.add_argument("--nuclear-fasta", required=True, type=Path, help="Nuclear assembly used for the catalog (its .fai gives the sequence lengths).")
    gi.add_argument("--gff", type=Path, help="Mitogenome annotation.")
    gi.add_argument("--long-read-bam", type=Path, help="Reads on the multiplied reference (redwood long_reads.redwood.bam).")
    gi.add_argument("--rnaseq-bam", type=Path)
    gi.add_argument("--variant-table", type=Path, help="redwood variants TSV for the variant ring.")
    go = parser_comp.add_argument_group("output")
    go.add_argument("--output-base", required=True, type=Path, help="Writes <base>.pdf/.png and <base>.legend.md.")
    go.add_argument("--fileform", nargs="+", default=["pdf", "png"])
    go.add_argument("--dpi", type=int, default=300)
    gs = parser_comp.add_argument_group("layout and style")
    gs.add_argument("--style", choices=sorted(STYLES), default="nature-communications")
    gs.add_argument("--layout", choices=["figure", "page"], default="page",
                    help="figure: just the figure at the figure width; page: the figure on a page with the legend underneath (default).")
    gs.add_argument("--page", help="Page size: letter, a4, nature-communications or WxH with a unit (210x279mm, 8.5x11in). Default: the style's page.")
    gs.add_argument("--figure-width", help="Figure width, e.g. 170mm or 6.7in. Default: the style's width (nature-communications: 170 mm).")
    gs.add_argument("--legend-columns", type=int, choices=[1, 2], help="Legend columns (nature-communications: 2).")
    gs.add_argument("--panel-labels", choices=["lower", "upper"], help="Panel letters a-d or A-D (default: the style's).")
    gs.add_argument("--panel-label-size", type=float, help="Panel letter size in pt (nature-communications: 8).")
    gs.add_argument("--text-size", type=float, nargs=2, metavar=("MIN", "MAX"), help="Figure text size range in pt (nature-communications: 5 7).")
    gs.add_argument("--font", action="append", help="Font family to try first (repeatable). The style's list follows (Helvetica, Arial, ...).")
    gs.add_argument("--font-dir", action="append", help="Directory of .ttf/.otf fonts to register (repeatable).")
    gl = parser_comp.add_argument_group("label and legend")
    gl.add_argument("--figure-label", help='Exact label, e.g. "Figure S4" or "Supplementary Fig. 4" (overrides the options below).')
    gl.add_argument("--figure-number", help="Figure number; the label is the style's prefix + number.")
    gl.add_argument("--supplementary", action="store_true", help="Use the style's supplementary prefix (nature-communications: Supplementary Fig.).")
    gl.add_argument("--label-prefix", help='Prefix word to use with --figure-number instead of the style\'s, e.g. "Supplementary Figure".')
    gl.add_argument("--species", help="Species name for the default title (set in italics).")
    gl.add_argument("--title", help="Title sentence after the label; markup: **bold**, *italic*.")
    gl.add_argument("--caption-file", type=Path, help="Legend body with markup and {fields}; default: a generated description of panels a-d.")
    gl.add_argument("--caption-append", help="Sentence(s) appended to the legend body.")
    gl.add_argument("--field", action="append", metavar="KEY=VALUE", help="Extra {KEY} value for the legend text (repeatable).")
    gl.add_argument("--no-legend", action="store_true", help="No legend on the page and no <base>.legend.md.")
    parser_comp.set_defaults(func=run_composite)

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
    parser_run.add_argument("--min-minor-frac", type=float, default=0.05,
                            help="Variant table / marks: flag columns whose minor-allele or insertion fraction reaches this.")
    parser_run.add_argument("--min-base-quality", type=int, default=20, help="Variant table: minimum base quality.")
    parser_run.add_argument("--read-mismatches", dest="read_mismatches", choices=["shared", "all", "none"], default="shared",
                            help="IGV-style mismatch marks on the read rings (see `redwood plot --help`).")
    parser_run.add_argument("--no-variant-ring", dest="no_variant_ring", action="store_true")
    parser_run.add_argument("--no-track-legend", dest="no_track_legend", action="store_true", help="Omit the ring key on the circular plot.")
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
