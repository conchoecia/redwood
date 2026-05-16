# redwood

`redwood` is a standalone circular genome plotting tool extracted from the
`pauvre redwood` plotter. It draws circular plots with optional long-read BAM
rings, GFF annotation tracks, and RNA-seq depth tracks.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/assets/redwood-grid-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="docs/assets/redwood-grid-light.png">
  <img alt="Four example redwood circular genome plots for human, mouse, fly, and sponge mitogenomes." src="docs/assets/redwood-grid-light.png">
</picture>

## Install

```bash
pip install .
```

For development:

```bash
pip install -e ".[dev]"
```

## Usage

Plot annotation only:

```bash
redwood plot --gff tests/testdata/gff_files/Bf201706.gff \
  --no-timestamp -o Bf201706
```

Plot long reads, annotation, and RNA-seq depth:

```bash
redwood plot \
  --main-bam reads.bam \
  --rnaseq-bam rnaseq.bam \
  --gff annotation.gff \
  --doubled main \
  --query "ALNLEN >= 10000" "MAPLEN < reflength" \
  -o sample_redwood
```

Input BAM files must be indexed with `samtools index`.

## End-to-End Workflow

`redwood run` builds the intermediate references and BAM files needed for a
complete circular genome plot. Provide a mitochondrial genome, optional GFF3
annotation, long reads, and RNA-seq reads. If RNA-seq reads are supplied, also
provide a nuclear genome FASTA; mitochondrial-looking nuclear FASTA contigs are
excluded before the mitochondrial FASTA is appended as the target.

```bash
redwood run \
  --mito-fasta mitochondrion.fa \
  --nuclear-fasta nuclear.fa \
  --gff annotation.gff \
  --long-reads ont.fastq.gz \
  --rnaseq-reads rna_1.fastq.gz rna_2.fastq.gz \
  --outdir redwood-work
```

The workflow maps long reads to a doubled mitochondrial reference, selects reads
that span large fractions of the circular genome, maps RNA-seq reads against a
nuclear bait reference plus the mitochondrial genome, keeps primary
mitochondrial RNA-seq alignments, and writes `redwood.metrics.json`,
`redwood.workflow.json`, and a plot output base under the requested output
directory. It expects `minimap2` and `samtools` on `PATH`.

For existing BAM files, use `redwood plot` directly. To summarize existing BAMs
without plotting, use `redwood metrics`.

The plotting CLI also accepts `--extra-track` declarations for newer plot
styles, including `at`, `gc`, `rnaseq-strand`, and `metrics`. The legacy plotter
currently renders the default read, annotation, and RNA-seq depth tracks.

## Notes

This repository keeps the original redwood plotting lineage from
[`pauvre`](https://github.com/conchoecia/pauvre), focused into a dedicated
package and command-line interface for circular genome plots.
