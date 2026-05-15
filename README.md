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

## Notes

This repository keeps the original redwood plotting lineage from
[`pauvre`](https://github.com/conchoecia/pauvre), focused into a dedicated
package and command-line interface for circular genome plots.
