# Example datasets

These examples are intended for visual work on redwood plots. BAM files in this
directory are committed fixtures from real sequencing runs, not simulated reads.
Most are trimmed to roughly 25-40x mitochondrial depth so plot iteration stays
fast, with reads selected for long mitochondrial spans and balanced circular
coverage.

Build or refresh the circular datasets with:

```bash
python scripts/build_real_read_datasets.py
```

The builder downloads a bounded prefix of each run's public FASTQ file(s), maps
those real reads to the target mitogenome, and keeps only mapped alignments.
Temporary FASTQ files stay under `examples/datasets/.cache/` and are not
committed. Long-read fixtures are selected greedily to favor reads that span a
large fraction of the mitogenome while filling under-covered circular bins.

## Beroe manuscript reads

Schultz et al. 2020, PeerJ (`PMC6991124`) states that Beroe mitochondrial reads
are available through BioProject `PRJNA421807`. NCBI's BioProject page currently
links BioSamples with SRA sample identifiers (`SRS2786396`, `SRS2786397`,
`SRS2786398`, `SRS2786399`, `SRS5502111`), but NCBI SRA runinfo, NCBI SRA
search, and ENA run search did not expose run accessions for those sample IDs
on 2026-05-15. The builder records this in `beroe/manifest.json` and does not
fabricate Beroe BAMs.

## Included real-run targets

- `human`: `NC_012920.1`, real ONT WGS reads from `DRR165688`.
- `drosophila`: `NC_024511.2`, real ONT WGS reads from `SRR12187559`.
- `mouse`: `NC_005089.1`, real ONT WGS reads from `DRR188136`.
- `sponge`: `NC_010202.1`, real PacBio WGS reads from `SRR10983242`
  (`Ephydatia muelleri`, complete annotated mitogenome).
- `hydra`: `NC_010214.1`, real ONT WGS reads from `SRR18362010`
  (`Hydra oligactis`, linear mitochondrial reference).

Each dataset directory contains a `manifest.json` with source run metadata and
a command to render the plot.

The circular BAMs are mapped to a doubled copy of each mitogenome so reads crossing the
linearized origin can be plotted correctly. Use `--doubled main` when rendering
datasets whose manifest has `"doubled_reference": true`.

## Linear Hydra example

The linear example uses the public *Hydra oligactis* mitochondrial reference
[NC_010214.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_010214.1), 16,314 bp, and ONT
run [SRR18362010](https://www.ncbi.nlm.nih.gov/sra/SRR18362010) from
[PRJNA816482](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA816482). The sequencing
study is [Cazet et al., Genome Research (2023)](https://doi.org/10.1101/gr.277040.122).
The published reference and WGS reads are from different source specimens.
This is a plotting example, not a newly polished assembly.

Rebuild it with `minimap2` on PATH and the ONT optional dependencies installed:

```bash
pip install -e '.[ont]'
python scripts/build_linear_read_dataset.py
```

The builder takes the first 10,000 records from the public FASTQ, maps them to the
single reference, extracts target-mapped reads, trims supported ONT end adapters
with Redwood's default Cutadapt settings, and remaps. The fixture retains all
mapped primary and supplementary records from that prefix. It does not perform
nuclear co-mapping or assign mitochondrial/NUMT classes. The README displays up
to 30 terminal-balanced reads; its depth, endpoint, and stacked soft-clip tracks
use all reads in the fixture, not the entire public sequencing run. No RNA-seq
track is shown for this dataset. The public reads commonly align to within
about 70–100 bp of the reference's right edge, so this example uses
`--terminal-window 100` to recognize both ends during selection. This changes
the end-group definition, not alignment coordinates or the reference sequence.

Gene/RNA coordinates come from the public RefSeq GFF. Display labels use gene
symbols (and amino-acid names for tRNAs); the partial `COX1_C` pseudogene is
retained as a gene arrow with its `pseudo=true` qualifier. Two additional
`repeat_region` features mark the **1,274 bp exact inverted-repeat core** at
68–1,341 and 14,827–16,100 (1-based, inclusive). The builder identifies this core
as the longest exact match between the first 2,000 bases and the reverse
complement of the last 2,000 bases. These are sequence-derived core boundaries,
not a claim about the full biological repeat extent. Terminal caps are not
annotated in this public example.

The [Hydra manifest](hydra/manifest.json) records source and fixture checksums,
mapping options, tool versions, read counts, and adapter-processing parameters.
[trimming.tsv](hydra/trimming.tsv) records removed sequences and lengths for
each candidate read. Downloaded prefixes and complete processing reports stay
in the ignored cache directory. Re-rendering the committed fixture with
`python scripts/build_readme_figure.py --only linear` needs neither downloads nor
Cutadapt/minimap2.

## RNA-seq bait mapping

The RNA-seq track should be built by mapping public RNA-seq reads against a
competitive reference containing the nuclear genome as bait plus the
mitochondrial genome as the target. That lets NuMT-like reads choose their best
nuclear placement instead of inflating mitochondrial expression.

`scripts/build_rnaseq_tracks.py` records one public RNA-seq candidate per
fixture species and expects a local whole-genome FASTA for the bait reference:

```bash
python scripts/build_rnaseq_tracks.py human --bait-fasta GRCh38.primary_assembly.fa
```

The script removes mitochondrial-looking contigs from the bait FASTA, appends
the fixture mitogenome, maps a bounded FASTQ prefix with `minimap2 -ax sr`, and
keeps only primary alignments whose best target is the mitochondrial contig.
For the committed RNA-seq BAMs, prefixes were increased until mean mitochondrial
depth was around 50x or higher; the fly RNA-seq BAM was then downsampled because
the initial prefix was far deeper than needed for plotting.

The current RNA-seq fixtures are:

- `human`: `DRR001175`, Illumina single-end RNA-seq.
- `mouse`: `DRR001494`, Illumina single-end RNA-seq.
- `drosophila`: `DRR016419`, Illumina single-end RNA-seq.
- `sponge`: `SRR14102585`, Illumina paired-end RNA-seq from `Ephydatia muelleri`.

If we use long-read RNA-seq later, pass `--preset splice` or another minimap2
preset. These candidates are short-read RNA-seq, so `minimap2 -ax sr` is the
appropriate starting point. For sponge intron/splicing questions, we should
compare this against a splice-aware short-read mapper before treating gaps as
biology.
