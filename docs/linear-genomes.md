# Linear mitogenome figures

Use `--topology linear` to draw a horizontal molecule with independent left and
right ends. A single-sequence `--mito-fasta` supplies the reference name, length,
and first/last 30 bases. BAM and GFF reference names must match it exactly, and
the BAM reference length must equal the FASTA length. A doubled BAM or
`--doubled` is rejected. Circular plotting remains the default.

## Figures from existing data

```bash
redwood plot --topology linear \
  --mito-fasta mitochondrion.fa \
  --gff annotation.gff3 \
  --main-bam long_reads.bam \
  --read-classes reads.tsv \
  --max-reads 80 \
  --title 'Linear mitochondrial genome' \
  --fileform pdf svg png --dpi 300 --no-timestamp -T \
  -o figures/mitochondrion
```

PDF and SVG preserve vector tracks, and SVG retains editable text. `-T` gives
an opaque background; the CLI otherwise defaults to transparency. `--dark`
selects a dark theme. `--width` controls figure width in inches (default 13).

For a compact annotation figure, omit `--main-bam` and `--read-classes`.
For annotation and depth without individual reads, use `--max-reads 0`.
`--hide-evidence` hides the endpoint and clip panels while keeping their exports.
An optional `--rnaseq-bam` adds primary RNA alignment depth; add
`--extra-track rnaseq-strand` to show the strands separately.

Gene/CDS/RNA arrows sit above or below the reference by strand. Redundant gene
parents at the same interval as a CDS or RNA are suppressed. Standalone genes,
including partial ITR copies, remain visible. `repeat_region` and `misc_feature`
annotations have separate lanes; inverted repeats with opposite GFF strands
have mirrored arrows. Features that overhang a terminus are clipped at the plot
boundary without wrapping to the other end. The original coordinates are not
changed. A GFF `##sequence-region`, when present, must describe the whole FASTA
sequence. Variant markers and haplotype panels are not part of this view.

## Read classes and displayed reads

The class file is a tab-separated table with a header:

```text
read	class	locus
read_001	mito	.
read_002	NUMT_flanked	chr2:1000-18000
read_003	chimera	chr7:8000-9000
```

`locus` is optional; `nuclear_locus` is also accepted, so the original
`phase_reads_fast.py` tables can be used directly. Additional columns are
ignored. Missing read IDs are `unclassified`. Contradictory duplicate labels
are rejected. Classes color the read panel and the per-class depth curves;
legends report **unique read counts**, not alignment counts. These are supplied
labels, not classifications inferred by Redwood. Recent NUMTs identical to
mtDNA cannot be separated by SNPs alone.

All input primary and supplementary alignments contribute to depth and endpoint
tracks. Secondary alignments are counted and excluded. A read contributes at
most once per reference base across overlapping segments. Depth counts aligned
bases (`M`, `=`, `X`), excluding deletions and reference skips; it has no pileup
depth cap. `--depth-scale log` uses a `log(1 + depth)` transform so zero-depth
positions remain visible.

`--max-reads` limits **only the read panel**, choosing the longest alignments by
default and grouping the selected reads by class. Supplementary segments of a
selected read share its row. This is a display subset, not a population-frequency
estimate. `--sort POS` sorts by increasing start instead. Insertions and
deletions use `--min-indel` (default 10 bp); soft clips are dots at the aligned
endpoints, without inventing reference sequence outside the molecule.

Explicit `--query` clauses retain the legacy column definitions:

| Column | Meaning |
| --- | --- |
| `POS` | 1-based alignment start |
| `MAPQ` | Mapping quality |
| `ALNLEN` | Sum of CIGAR lengths except `S`, `H`, `I` |
| `MAPLEN` | Sum of CIGAR lengths except `I` |
| `TRULEN` | Sum of all CIGAR lengths (legacy definition, including deletions) |
| `reflength` | Single-copy reference length |
| `READ`, `CLASS` | Read name and supplied class |

For example, `--query "ALNLEN >= 10000" "MAPLEN <= reflength"` applies both
clauses to the displayed alignment rows. Depth and evidence are unaffected by
queries and display limits. Linear plots default to no query, so full-length
reads and short reads on small genomes are retained. Use a BAM subset that
retains supplementary alignments if the evidence should describe one population.

## Terminus evidence and exports

For every linear plot with a main BAM, the same output prefix receives:

* `.evidence.json`: alignment/read counts, per-class read counts and mean depth,
  terminal endpoints, clip counts, raw junction counts, read-length median/p99,
  fraction with at least 300 bp unaligned to mitochondrial segments, and
  terminal 2% / middle 50% mean-depth ratios.
* `.profiles.tsv`: unbinned per-base class depths, alignment starts/ends,
  short/long clip counts on each side, and MAPQ-0 start/end positions.
* `.junctions.tsv`: adjacent segments ordered along the original read, their
  reference coordinates, strand, read class, locus, query gap, and join category.

TSV positions are 1-based; junction endpoints are inclusive. Endpoint tracks
count **alignments**, including supplementary segments. The uniform endpoint
expectation is therefore `alignment_count * window / reference_length`.
They show reference-coordinate left/right endpoints, independent of read strand.
Plot bars aggregate `--bin-size` bases (default 25) on a symmetric logarithmic
axis; exact terminal counts use `--terminal-window` (default 30 bp).

Clips of 1–99 bp and at least 100 bp are separate categories by default;
`--clip-threshold` changes that cutoff. These are length categories only.
Redwood does not establish cap identity or map clipped sequence to nuclear DNA.
Hard-clipped terminal alignments are marked unassessed in the summary.

Split joins use query-facing breakpoints within `--junction-window` (default
300 bp) of a terminus. Categories are end-to-start on the same strand, terminal
inverted, internal same-strand, and internal inverted. Gaps or overlaps exceeding
500 query bases and hard-clipped pairs are explicitly unassessed. Read orientation
is accounted for on both strands. These are **candidate alignment geometries**,
not diagnoses of circles, hairpins, deletions, or recombination.

Warnings cover hard clips, absent supplementary records, soft masking, and a
high unaligned-read fraction. With no supplementary records, zero split joins
is not evidence that joins were absent before BAM filtering. High unaligned
fractions can reflect NUMTs or chimeras and need nuclear co-mapping to interpret.

`--topology linear` selects the drawing and mapping model. It does **not** infer
topology; `topology_call` is `not_assessed`. The full `redwood topology` scorecard,
rotated-reference mapping, self-alignment, clip-content classification, and
automatic biological verdicts in [issue #3](https://github.com/conchoecia/redwood/issues/3)
remain separate work.

## Mapping workflow

```bash
redwood run --topology linear \
  --mito-fasta mitochondrion.fa --gff annotation.gff3 \
  --long-reads ont.fastq.gz --long-read-preset map-ont \
  --outdir linear-work --fileform pdf png
```

Linear mapping uses one copy of the reference and
`minimap2 -a -x PRESET -Y --secondary=no`. All primary and supplementary
alignments are retained. Circular read-selection thresholds (`--long-read-depth`
and `--min-span-fraction`) do not apply; use `--max-reads` to control the figure.
No circular junction or tandem reference is generated.

The same topology option is available to `advanced prepare-reference`,
`advanced map-long`, and `metrics`. Linear `map-long` rejects `--copies` greater
than one. The existing RNA bait workflow is unchanged. Long-read nuclear
co-mapping and automatic NUMT inference are not implemented in this PR; use
existing class tables for population-aware figures.
