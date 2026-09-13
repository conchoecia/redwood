# Working with redwood on a new mitogenome (for agents and people)

This file is the operating manual for using redwood to decide whether a candidate mitochondrial genome is right and whether
the reads that disagree with it are mitochondrial (heteroplasmy, a second lineage in the sample) or nuclear (NUMTs).
It distills the Hydra oligactis and clitellate (Helobdella, Haemadipsa, Perionyx) mitogenome tracks.

## 0. Ground rules
- An assembler's contig is a hypothesis. Never call a mitogenome "circular", "complete" or "done" until the reads say so.
- Bait reads against **mitochondrial sequence only**. A nuclear flank in the bait pulls in tens of thousands of repeat reads.
- Keep the whole read sequence in every alignment (`minimap2 -Y --secondary=no`) and keep supplementary alignments
  (`samtools view -F 0x104`); classification needs to see where the *rest* of a read went.
- Long reads decide structure (duplications, repeat copy number, circularity); accurate reads (HiFi / Illumina) decide bases.
- Write down which reads support each claim, and hand over BAMs. Nobody should have to trust a summary.

## 1. Candidate molecule
`MitoHiFi -r` (HiFi), Flye (CLR/ONT) or a reference-guided consensus give a candidate. Rotate it to a conventional start (COX1).
Expect these failure modes and test for them in step 2: a collapsed tandem duplication (2x read depth over a block, split reads
joining the block to itself, several assembler contigs carrying the block twice), a mosaic control region (assembler stitches
repeat arrays of different molecules), a NUMT contig instead of the molecule (a contig that runs into a nuclear flank; depth drops
from thousands to tens inside the contig), a linear molecule forced into a false circle through a terminal repeat.

## 2. Circularity and structure from the reads
```
redwood run --mito-fasta cand.fa --long-reads reads.fq.gz --long-read-preset map-hifi --outdir rw   # doubled/tripled reference
redwood metrics ... ; redwood variants --mito-fasta cand.fa --bam rw/long_reads.raw.bam --output variants.tsv
```
Read the metrics: `origin_spanning_alignments` at the depth of the rest of the molecule and a flat depth through the junction
mean a circle; a spike of alignment starts/ends at one internal coordinate means a terminus or a mis-join; a region at 2x depth
means a collapsed duplication; `redwood variants` columns tagged `mismatch` mean the candidate is not the read majority.
For a linear molecule use `--topology linear` (terminus pileups, soft-clip caps, split-read joins).
Measure structural features **inside reads** rather than trusting alignments: the distance between two unique anchors flanking a
repeat array / duplication, measured per read, gives the population's length distribution directly.

## 3. NUMTs: catalog, classify, plot
```
redwood numts --mito-fasta final.fa --nuclear-fasta nuclear.fa --gff final.gff --long-reads reads.fq.gz \
              --mito-bam rw/long_reads.raw.bam --outdir numts --figure --figure-bam rw/long_reads.redwood.bam
```
Outputs: `numts.loci.tsv` (locus, chromosome, coordinates, span, mitogenome intervals, fraction of the molecule, identity, class),
`numts.mito_coverage.tsv`, `numts.read_classes.tsv`, `numts.junction_support.tsv`, `reads_vs_mito.PO.bam` (reads tagged `PO:Z`),
`numts.landscape.png` (chromosomes with loci at their true span, widened to a stated minimum when too narrow to see,
colored by identity), `numts.catalog.png` (loci over the mitogenome at their
identity, alignments of one locus joined by a dotted line), `numts.sizes.png` (mtDNA content per locus against the class
thresholds: `full-length` >= 95 % of the molecule, `large` >= 5 kb, otherwise `fragment`), and with `--figure` the composite
one-page figure, 6.5 x 9 in to fit the text area of letter paper (a circular map with its ring key, b class panel, c landscape,
d catalog). Every circular map carries the ring key in its bottom-right corner (`--no-track-legend` to omit). `redwood plot --numt-loci numts.loci.tsv
--circular-read-classes numts.read_classes.tsv` adds the NUMT ring and colors reads by class on the circular map.

**Before cataloging, remove mitochondrial scaffolds from the nuclear assembly** (a standalone scaffold that is a jumble of
mitogenome copies is the assembler's mito output, not a NUMT; leaving it in makes every mito read "nuclear_only_at_NUMT_locus").

Read classes and how to read them. They describe where the segments of a read align, not proven origin:
| class | meaning |
|---|---|
| `mito_only` | one alignment, entirely on the mitogenome |
| `mito_multisegment` | two or more mito segments and nothing else: a read crossing the origin of a circular molecule (or a concatemer / rearrangement; check the join) |
| `mito+nuclear_at_NUMT_locus` | candidate junction read: a nuclear segment next to a cataloged NUMT with >= 500 bp (`--min-flank`) outside the locus; confirm the query order and the junction before calling the read NUMT-derived |
| `mito+NUMT_homology` | every nuclear segment lies mostly inside cataloged NUMT sequence: a competing placement of mitochondrial-like sequence, no nuclear flank; usually a mitochondrial read |
| `mito+nuclear_elsewhere` | mito + nuclear segments away from cataloged NUMTs: an uncataloged or assembly-absent NUMT if independent reads share the nuclear breakpoint, otherwise a library chimera |
| `nuclear_only_at_NUMT_locus` | every segment nuclear and overlapping a cataloged locus |
| `nuclear_only` | nuclear segments only; mostly bait carry-over with a short mito-like hit, but an assembly-absent NUMT also lands here when its mitochondrial sequence sits in an unaligned insertion or clip (check recurrent insertions/clips at one breakpoint before discarding) |

## 4. Is a divergent read population mitochondrial or NUMT? Decision rules
1. **Flanks first.** Reads with a substantial, uniquely placed nuclear flank at a locus, joined to the mitochondrial part in a
   compatible query order, are NUMT-derived; count them per locus and confirm the insertion with `numts.junction_support.tsv`
   (reads with aligned bases on both sides of each boundary). A "flank" shorter than ~500 bp, or one that lies inside a short NUMT
   fragment 98-99 % identical to mtDNA (`mito+NUMT_homology`), is a competing placement, not a flank.
2. **Reads longer than any NUMT.** If the longest cataloged NUMT is 8 kb and the divergent reads are 15-25 kb of contiguous
   mitogenome, they are not explained by a cataloged NUMT. The catalog only covers what the assembly contains and the search
   thresholds detect: look for uncataloged or assembly-absent NUMTs via `mito+nuclear_elsewhere` reads, and recurrent
   insertions or clips carrying mitochondrial sequence, that share one nuclear breakpoint.
3. **Origin-crossing reads support, but do not prove, mitochondrial origin.** The origin of the reference is arbitrary, so a
   single nuclear insertion can contain it, and a tandem nuclear copy need not leave a spacer at the join. Many divergent reads
   that wrap the circle with long continuous alignments, compatible joins and linked alleles, and no nuclear flanks, are best
   explained as mitochondrial (a second lineage or heteroplasmy); state the conclusion as "not explained by the cataloged
   nuclear loci".
4. **Repeat copy-number variation fits a mitochondrial origin.** Indels in reads that are integer multiples of the control-region
   repeat unit and vary between reads fit heteroplasmic VNTR length variation, but NUMTs can also differ between loci, alleles or
   pooled individuals, so use it together with rules 1-3, not alone.
5. **SNV phasing separates populations only with accurate reads.** With HiFi, cluster reads by their alleles at variant sites
   (minor allele >= a few %): mitochondrial lineages appear as clusters differing at many linked sites; NUMT reads form clusters
   that coincide with flanked reads. With CLR/ONT (error > NUMT divergence) rely on rules 1-4 and on `--read-mismatches shared`
   with the error-adaptive threshold.
6. **Restrict phasing sites to unique sequence.** Sites inside repeat arrays / duplications are alignment noise; if most "variant
   sites" fall in a repeat block, reads that do not cover the block cannot be assigned and will look like a haplotype gap.
7. **Fresh NUMTs identical to mtDNA cannot be separated by sequence**, only by flanks and depth. Say so.
   Report NUMT counts as "detected in this assembly under these thresholds" (minimap2 asm20 + dc-megablast >= 150 bp, 3-kb
   merge), and give locus span and mtDNA content separately: a merged locus is not necessarily one insertion event.
8. **Where the mitogenome disagrees with its own reads** (`mismatch` columns in `redwood variants` computed on the population
   you believe is the molecule), rebuild the consensus from that population only (`samtools consensus -c 0.51` on their BAM);
   racon on HiFi collapses repeat arrays and on CLR shrinks the molecule, so use it only to polish a de novo CLR/ONT contig.

**Figures for a paper:** `redwood composite --style nature-communications` (or `generic`) lays out the map and NUMT panels to the journal's measured geometry, as a figure (`--layout figure`) or a full page with the legend (`--layout page`); set `--page`, `--figure-width`, `--panel-labels`, `--figure-label` or `--figure-number --supplementary` as the manuscript needs, and edit `<base>.legend.md` rather than retyping the numbers.

## 5. Deliverables that let someone else check the claim
Final FASTA (rotated to COX1) + GFF/GenBank; the doubled/tripled-reference BAM and a single-copy BAM with `HP` (haplotype cluster)
and `PO` (segment class) tags; per-population BAMs; `variants.tsv`; `numts.loci.tsv` + `numts.read_classes.tsv` +
`numts.junction_support.tsv`; the redwood map, the NUMT landscape and catalog; an IGV session that opens them; a README with
the counts behind every statement (origin-crossing reads, junction depth ratio, reads per lineage, loci per class).
