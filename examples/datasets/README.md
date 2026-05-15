# Example datasets

These examples are intended for visual work on redwood plots. BAM files in this
directory are committed fixtures from real sequencing runs, not simulated reads.
Most are trimmed to roughly 30-50x mitochondrial depth so plot iteration stays
fast.

Build or refresh datasets with:

```bash
python scripts/build_real_read_datasets.py
```

The builder downloads a bounded prefix of each run's public FASTQ file(s), maps
those real reads to the target mitogenome, and keeps only mapped alignments.
Temporary FASTQ files stay under `examples/datasets/.cache/` and are not
committed.

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
- `drosophila`: `NC_024511.2`, real ONT WGS reads from `SRR12177581`.
- `mouse`: `NC_005089.1`, real ONT WGS reads from `DRR188136`.
- `sponge`: `MZ675556.1`, real Illumina WGS reads from `SRR15070719`
  (`Spheciospongia vesparium`, 21.8 kb mitogenome).

Each dataset directory contains a `manifest.json` with source run metadata and
a command to render the plot.

The BAMs are mapped to a doubled copy of each mitogenome so reads crossing the
linearized origin can be plotted correctly. Use `--doubled main` when rendering
datasets whose manifest has `"doubled_reference": true`. The sponge dataset is
about 26x because the mitochondrial fraction in the bounded read prefix is low.
