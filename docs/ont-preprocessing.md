# ONT adapter trimming and remapping

Redwood can reprocess the complete ONT molecules that have an alignment to the
target genome. It extracts those reads, searches their ends with Cutadapt,
trims supported adapters, and remaps all extracted molecules with minimap2
`map-ont`. Plotting and evidence then use the new indexed BAM.

Install the optional Python dependency and put minimap2 on `PATH`:

```bash
pip install 'redwood[ont]'
# From a checkout:
pip install -e '.[ont]'
```

The mapping workflows also require samtools. Plain plotting of an existing BAM
does not require Cutadapt or minimap2 unless reprocessing is requested.

## Existing ONT BAM

Declare that the main BAM contains ONT reads when plotting:

```bash
redwood plot --topology linear --reprocess-ont \
  --mito-fasta mitochondrion.fa --main-bam ont.bam --gff genes.gff3 \
  --linear-track depth --linear-track ends --linear-track clips \
  --terminal-details --fileform pdf svg png --no-timestamp -T \
  -o figures/ont_trimmed
```

This writes the preprocessing files to `figures/ont_trimmed.ont-preprocess/`.
Choose another **new** directory with `--ont-workdir`. Redwood refuses to
overwrite an existing preprocessing directory. To redraw without processing
again, use its `reads.trimmed.remapped.bam` as `--main-bam` and omit
`--reprocess-ont`; retain the preprocessing report with the figure's methods.

For preprocessing without a figure:

```bash
redwood advanced reprocess-ont --topology linear \
  --mito-fasta mitochondrion.fa --main-bam ont.bam --outdir ont_processed
```

Redwood includes each read ID once if any of its alignments maps to the target,
including supplementary or secondary alignments. It restores original read
orientation and keeps soft-clipped DNA. Full sequences and qualities must be
recoverable from the BAM; for hard-clipped or sequence-missing records, provide
`--ont-source-reads original.fastq.gz`. Conflicting full records with the same
read ID cause an error. The target FASTA must contain one sequence. Linear
mode requires matching BAM reference name and length. Circular mode preserves
the existing number of tandem reference copies.

This is candidate extraction, not mitochondrial-versus-NUMT classification.
Reference-only remapping cannot establish that a candidate is mitochondrial.
Existing read-class tables remain labels keyed by read ID; Redwood does not
recalculate their classifications after trimming.

## Workflows starting from reads

`redwood run` and `redwood advanced map-long` automatically perform this step
after initial mapping when their long-read preset is `map-ont`. Other presets,
including `map-hifi` and `map-pb`, do not trigger ONT trimming. Use
`--no-ont-trim` to skip it for already processed data. ONT input must be FASTQ
with base qualities. For circular genomes, trimming/remapping precedes the
existing circular read selection. Linear mode keeps all remapped candidates
and limits only the displayed read sample.

## Adapter matching

Defaults are two published ONT ligation adapter motifs and their corresponding
reverse complements, from [Porechop's adapter definitions](https://github.com/rrwick/Porechop/blob/master/porechop/adapters.py).
This default does not identify a sequencing kit or cover every ONT chemistry.
Provide library-specific sequences with repeated `--ont-adapter-5p` and
`--ont-adapter-3p` options. Each option replaces the defaults for that end;
sequences must be A/C/G/T in the original FASTQ orientation.

- `--ont-overlap 12`: minimum matching adapter overlap in bases.
- `--ont-error-rate 0.1`: maximum substitution/indel error rate.
- `--ont-end-window 100`: search only the outer 100 bases at each end, capped
  at half the read length so the two search windows do not overlap.
- `--ont-threads 4`: trimming and remapping threads.

Redwood passes the two end windows separately to Cutadapt, allowing partial
adapters at the actual read ends and full matches with flanking noise inside
the window. A match removes the adapter plus sequence toward the outside read
end. These boundaries are applied to the complete original molecule, preserving
its orientation, ID and remaining qualities. See the
[Cutadapt adapter matching guide](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types).
These conservative defaults are a starting point to inspect, not a validated
optimum for every library. Adapters beyond the end window are left for separate
chimera analysis. No quality trimming, G/C-run trimming, read-length filtering,
or automatic splitting of internal adapters occurs. A candidate that would
become empty is retained unchanged and flagged in the trimming table.

## Audit outputs

- `reads.original.fastq.gz`, `reads.trimmed.fastq.gz`: complete candidate reads
  before and after trimming, including candidates that fail to remap.
- `removed_tails.fastq.gz`, `trimming.tsv`: exact removed sequences, lengths,
  original qualities and read IDs, in original sequencing orientation. The
  table also counts removed bases previously aligned in the primary target
  record; a blank count indicates an incompatible original read length.
- `5p.cutadapt.json`, `3p.cutadapt.json`, corresponding logs/info TSVs and end
  FASTQs: adapter matches and tool output for both ends.
- `mapping_reference.fa`, `reads.trimmed.remapped.bam` and its index: the
  mapping reference and new primary/supplementary alignments. Unmapped
  candidates are retained in the BAM.
- `alignment_comparison.tsv`: per-read mapping status, primary alignment
  coordinates (1-based inclusive), soft clips and aligned-base counts before
  and after. Clipping columns use reference-left/right orientation.
- `ont_preprocessing.json`, `commands.json`, `minimap2.log`: parameters,
  tool versions, reference sequence checksum, counts and provenance.

Before/after differences include the effect of remapping; they should not all
be attributed to trimming if the original mapper/version/settings differed.
Inspect the removed sequences and retained terminal coverage before choosing
production inputs. Length-based display selection runs again after remapping
and may select different read IDs. When preprocessing runs as part of linear
plotting or `redwood run`, its report is also recorded in figure evidence and
the publication caption.
