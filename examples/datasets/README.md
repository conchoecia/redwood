# Example datasets

These examples are intended for visual work on redwood plots. BAM files in this
directory are committed fixtures from real sequencing runs, not simulated reads.
Most are trimmed to roughly 25-40x mitochondrial depth so plot iteration stays
fast, with reads selected for long mitochondrial spans and balanced circular
coverage.

Build or refresh datasets with:

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

Each dataset directory contains a `manifest.json` with source run metadata and
a command to render the plot.

The BAMs are mapped to a doubled copy of each mitogenome so reads crossing the
linearized origin can be plotted correctly. Use `--doubled main` when rendering
datasets whose manifest has `"doubled_reference": true`.

## RNA-seq bait mapping

The RNA-seq track should be built by mapping public RNA-seq reads against a
competitive reference containing the nuclear genome as bait plus the
mitochondrial genome as the target. That lets NuMT-like reads choose their best
nuclear placement instead of inflating mitochondrial expression.

`scripts/build_rnaseq_tracks.py` records one small public RNA-seq candidate per
fixture species and expects a local whole-genome FASTA for the bait reference:

```bash
python scripts/build_rnaseq_tracks.py human --bait-fasta GRCh38.primary_assembly.fa
```

The script removes mitochondrial-looking contigs from the bait FASTA, appends
the fixture mitogenome, maps a bounded FASTQ prefix with `minimap2 -ax sr`, and
keeps only primary alignments whose best target is the mitochondrial contig.
The initial candidates are:

- `human`: `DRR001175`, Illumina single-end RNA-seq.
- `mouse`: `DRR001494`, Illumina single-end RNA-seq.
- `drosophila`: `DRR016419`, Illumina single-end RNA-seq.
- `sponge`: `SRR3168560`, Illumina single-end RNA-seq from `Ephydatia muelleri`.

If we use long-read RNA-seq later, pass `--preset splice` or another minimap2
preset. These candidates are short-read RNA-seq, so `minimap2 -ax sr` is the
appropriate starting point. For sponge intron/splicing questions, we should
compare this against a splice-aware short-read mapper before treating gaps as
biology.
