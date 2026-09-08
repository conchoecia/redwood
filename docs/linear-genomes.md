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
`--extra-track rnaseq-strand` for strand coloring.

### Production style (default)

`--linear-style redwood` unrolls the existing circular plot's visual encodings:

* Bark-colored RNA depth (when supplied), with the same `log(1 + depth)` scaling.
* Green CDS/gene arrows, red rRNAs, pink tRNAs, and white labels inside features
  where they fit. Forward/reverse strands have separate lanes.
* The warm AT strip, using the circular renderer's colormap and robust 2nd–98th
  percentile color range. The 201 bp sequence windows shorten at linear termini
  instead of wrapping. `--extra-track gc` adds a GC strip.
* Tightly packed read rows with the original heartwood-to-sapwood gradient and
  the same insertion/deletion width ratios as regular circular reads. Each SVG
  read is one outline with a native continuous gradient; PDF uses vector shading
  clipped to that outline, instead of hundreds of solid color slices.
* The existing light/dark backgrounds, subtle frame, coordinate ticks, and
  centered reference-length label.

The production figure also includes **long-read depth** when a main BAM is
provided, and a thin **ITR/cap lane** from GFF repeat/misc-feature annotations.
The depth band shows total unique-read depth; per-class values are exported.
RNA and long-read depth are separate tracks with their own ranges.
Use `--rnaseq-label 'Sample RNA'` to identify the RNA sample in the production
figure. RNA depth counts aligned bases in primary records; barcode/UMI tags do
not trigger deduplication.

Add the full set of compact diagnostics in Redwood's bark colors:

```bash
redwood plot --topology linear --mito-fasta mitochondrion.fa \
  --gff annotation.gff3 --main-bam long_reads.bam \
  --linear-track depth --linear-track ends --linear-track clips \
  --fileform pdf svg png --no-timestamp -T -o figures/linear_evidence
```

`--linear-track` replaces the default depth-only selection; repeat it to select
several bands, or use `--linear-track none` for only the original core tracks and
terminal annotations. Start/end and clip counts use mirrored `log(1 + count)`
histograms, with scales and clip categories printed on the figure.
The short labels "Log depth" and "log count" refer to these `log(1 + value)`
transforms, which retain zeros. Printed maxima are raw depth/count values,
not transformed values; the endpoint maximum applies to either side of the
baseline and is a count per bin.
`--show-terminal-sequences` adds the first/last 30 bases below the plot.

### Diagnostic style

`--linear-style diagnostic` retains the separate annotation, read-class depth,
endpoint, and clip panels with explicit axes and full terminal sequences.
This view colors reads by the supplied classes. RNA strand mode shows separate
curves here. Both styles select the same reads and pack them by reference
position. The evidence calculations and exports are identical in both styles.

### Annotation handling

Gene/CDS/RNA arrows are separated by strand. Redundant gene
parents at the same interval as a CDS or RNA are suppressed. Standalone genes,
including partial ITR copies, remain visible. `repeat_region` and `misc_feature`
annotations are shown in terminal lanes; inverted repeats with opposite GFF strands
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
are rejected. In production figures, `--read-color class` replaces the wood
gradient with class colors and adds a count legend. Diagnostic figures color
reads and per-class depth curves automatically. Legends report **unique read
counts**, not alignment counts. These are supplied
labels, not classifications inferred by Redwood. Recent NUMTs identical to
mtDNA cannot be separated by SNPs alone.

All input primary and supplementary alignments contribute to depth and endpoint
tracks. Secondary alignments are counted and excluded. A read contributes at
most once per reference base across overlapping segments. Depth counts aligned
bases (`M`, `=`, `X`), excluding deletions and reference skips; it has no pileup
depth cap. `--depth-scale log` uses a `log(1 + depth)` transform so zero-depth
positions remain visible.

`--max-reads` limits **only the read panel**. Linear selection happens before
coordinate sorting, so taking the first displayed rows cannot bias the sample
toward the left terminus. The default `--linear-read-selection terminal-balanced`
selects the longest alignments in turn from three exclusive groups: left-only,
right-only, and spanning both termini. A read reaches a terminus when its
alignment boundary falls within `--terminal-window` bases (default 30).
Strand does not change which reference end is reached. Empty groups give their
slots to the remaining terminal groups; internal reads fill any slots left after
all terminal candidates are exhausted. A read is counted once, using its longest
primary alignment (or longest supplementary alignment if no primary is present).
Separate split segments do not qualify as a single end-to-end alignment.

`--linear-read-selection longest` selects the longest reads overall instead.
By default, length ranking uses `ALNLEN`; `--sort MAPLEN` or `--sort TRULEN`
changes the length measure used for selection. In linear mode, `--sort POS`
retains `ALNLEN` ranking rather than selecting a prefix of the coordinates.

After selection, both linear styles arrange reads by increasing leftmost
alignment start, then rightmost end, with read name breaking exact ties.
Nonoverlapping reads can share a row. All supplementary segments of a selected
read are retained on its row, even if a segment did not itself pass the query.
This is a display subset, not a population-frequency estimate. The evidence
JSON records the selection method, terminal groups, and selected read IDs.
Circular selection and ordering are unchanged.

Insertions and deletions use `--min-indel` (default 10 bp). Production reads encode these as
width changes, matching the circular plot. Diagnostic reads use insertion ticks,
deletion gaps, and soft-clip endpoint dots. Neither invents reference sequence
outside the molecule; the production clip band summarizes unaligned tails.
The reference runs left to right, from base 1 to the final base. The wood
gradient runs dark to light in that direction on both alignment strands; it is
not a read 5-prime/3-prime indicator. Display order and strand are independent.

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
clauses to read eligibility; a read qualifies when any segment passes all clauses.
Depth and evidence are unaffected by
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
These are alignment boundaries, not necessarily physical molecule ends. A
shorter read can end internally without clipping; a clipped or supplementary
alignment can end before the original read ends. Endpoint tracks include all
input alignments even when the read panel displays only the longest reads.
Plot bars aggregate `--bin-size` bases (default 25), using mirrored `log(1 + count)`
heights in production or symmetric logarithmic axes in diagnostic panels.
Exact terminal counts use `--terminal-window` (default 30 bp).

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
