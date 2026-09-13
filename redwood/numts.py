"""NUMT catalog, read classification and figures.

A nuclear genome usually carries copies of mitochondrial sequence (NUMTs). They confuse mitogenome assembly (an assembler may
extend the mitochondrion into a NUMT flank, or build the NUMT instead of the molecule) and read-population analyses (NUMT-derived
reads look like a divergent mitochondrial haplotype). This module

* catalogs NUMT loci: the mitogenome is aligned to the nuclear assembly (minimap2 ``asm20``; optionally ``blastn -task
  dc-megablast`` for short diverged fragments), hits on one chromosome within ``merge`` bp are merged into a locus, and each locus
  gets its span, the mitogenome intervals it covers, the fraction of the molecule covered, SNV-only and total identity, and a class
  (``full-length`` >= 95 % of the molecule, ``large`` >= 5 kb, ``fragment``);
* classifies long reads mapped to nuclear + mitochondrial sequence together (``minimap2 -Y --secondary=no``, supplementary
  alignments kept) by the origin of their segments: ``mito_only``, ``mito_multisegment`` (two mito segments, e.g. a read crossing
  the origin of a circular molecule), ``mito+nuclear_at_NUMT_locus`` (a nuclear flank of >= 500 bp at a cataloged NUMT),
  ``mito+NUMT_homology`` (nuclear segments only inside NUMT sequence: a competing placement, no flank), ``mito+nuclear_elsewhere``
  (an uncataloged NUMT or a library chimera), ``nuclear_only_at_NUMT_locus`` (a NUMT read whose whole length is explained by the
  nuclear locus) and ``nuclear_only``;
* counts reads spanning each nuclear-mitochondrial junction in a single alignment (junction support);
* draws the NUMT landscape (chromosomes with loci colored by identity), the NUMT-versus-mitogenome catalog (which part of the
  molecule each locus covers, at what identity), the per-locus mtDNA-content panel that defines the classes, and a composite
  letter-proportioned figure (circular map and class panel on top, landscape and catalog across the full width).
"""

from __future__ import annotations

import collections
import json
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam
from matplotlib.patches import Rectangle

CLASS_COLORS = {
    "mito_only": "#9c5a2c", "mito_multisegment": "#3a1d10", "mito+nuclear_at_NUMT_locus": "#d62728",
    "mito+nuclear_elsewhere": "#800080", "nuclear_only_at_NUMT_locus": "#7f7f7f", "nuclear_only": "#bdbdbd",
    "mito+NUMT_homology": "#c49a6c",
}
LOCUS_COLUMNS = ["locus", "chrom", "start", "end", "span_bp", "n_hits", "mito_bp", "nuclear_bp", "mtdna_bp", "frac_of_mito", "identity", "identity_snv",
                 "n_snv", "n_indel", "class", "strands", "mito_intervals"]


# ----------------------------------------------------------------------------------------------------------- catalog
def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else: buf.append(line.strip())
    if name is not None: seqs[name] = "".join(buf)
    return seqs


def _cs_counts(cs: str) -> tuple[int, int, int]:
    """(matched bp, substitutions, indel events) from a minimap2 cs string."""
    match = sum(int(x) for x in re.findall(r":(\d+)", cs))
    snv = len(re.findall(r"\*[a-z]{2}", cs))
    indel = len(re.findall(r"[+-][a-z]+", cs))
    return match, snv, indel


def minimap2_hits(mito_fasta: Path, nuclear_fasta: Path, threads: int = 8) -> list[dict]:
    """Mitogenome (query) vs nuclear assembly (target), all hits."""
    cmd = ["minimap2", "-t", str(threads), "-cx", "asm20", "--cs", "-N", "500", "-p", "0.01", "--secondary=yes",
           str(nuclear_fasta), str(mito_fasta)]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    hits = []
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 12: continue
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), "")
        nm = next((int(x[5:]) for x in f[12:] if x.startswith("NM:i:")), None)
        alen = int(f[10]); match, snv, indel = _cs_counts(cs)
        hits.append({"chrom": f[5], "tstart": int(f[7]), "tend": int(f[8]), "qstart": int(f[2]), "qend": int(f[3]),
                     "strand": f[4], "alen": alen, "identity": (alen - nm) / alen if nm is not None else int(f[9]) / alen,
                     "identity_snv": match / (match + snv) if match + snv else float("nan"), "snv": snv, "indel": indel, "source": "minimap2"})
    return hits


def blastn_hits(mito_fasta: Path, nuclear_fasta: Path, workdir: Path, threads: int = 8, min_len: int = 150,
                max_targets: int = 1_000_000) -> list[dict]:
    """Short / diverged fragments that asm20 chaining misses (dc-megablast). Needs BLAST+ on PATH; raises when it is
    missing, because a silently empty search would read as "no NUMTs"."""
    if not (shutil.which("blastn") and shutil.which("makeblastdb")):
        raise RuntimeError("BLAST search requested (--blastn) but blastn/makeblastdb are not on PATH")
    db = workdir / "nuclear_db"
    subprocess.run(["makeblastdb", "-in", str(nuclear_fasta), "-dbtype", "nucl", "-out", str(db)], check=True, capture_output=True)
    out = subprocess.run(["blastn", "-task", "dc-megablast", "-query", str(mito_fasta), "-db", str(db), "-evalue", "1e-10",
                          "-num_threads", str(threads), "-max_target_seqs", str(max_targets), "-max_hsps", "5000",
                          "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident"],
                         capture_output=True, text=True, check=True).stdout
    hits = []
    for line in out.splitlines():
        q, t, pid, ln, mm, go, qs, qe, ss, se, ev, bs, nid = line.split("\t")
        if int(ln) < min_len: continue
        ss, se = int(ss), int(se); strand = "+" if se >= ss else "-"
        hits.append({"chrom": t, "tstart": min(ss, se) - 1, "tend": max(ss, se), "qstart": int(qs) - 1, "qend": int(qe), "strand": strand,
                     "alen": int(ln), "identity": float(pid) / 100,
                     "identity_snv": int(nid) / (int(nid) + int(mm)) if int(nid) + int(mm) else float("nan"), "snv": int(mm),
                     "indel": int(go), "source": "blastn"})
    return hits


def _union(intervals: list[tuple[int, int]]) -> list[list[int]]:
    out: list[list[int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out


def merge_loci(hits: list[dict], mito_length: int, merge: int = 3000, full_length: float = 0.95, large_bp: int = 5000) -> list[dict]:
    """Cluster hits per chromosome (any hit within ``merge`` bp of the running locus joins it) and summarize each locus."""
    by_chrom: dict[str, list[dict]] = collections.defaultdict(list)
    for h in hits: by_chrom[h["chrom"]].append(h)
    loci = []
    for chrom, hs in by_chrom.items():
        hs.sort(key=lambda h: h["tstart"]); cur = [hs[0]]
        for h in hs[1:]:
            if h["tstart"] <= max(x["tend"] for x in cur) + merge: cur.append(h)
            else: loci.append((chrom, cur)); cur = [h]
        loci.append((chrom, cur))
    rows = []
    for chrom, hs in loci:
        start, end = min(h["tstart"] for h in hs), max(h["tend"] for h in hs)
        iv = _union([(h["qstart"], h["qend"]) for h in hs]); cov = sum(e - s for s, e in iv)
        weighted = [(h["identity"], h["alen"]) for h in hs if h["source"] == "minimap2" and h["alen"] >= 500] or [(h["identity"], h["alen"]) for h in hs]
        ident = sum(i * w for i, w in weighted) / sum(w for _, w in weighted)
        snvw = [(h["identity_snv"], h["alen"]) for h in hs if h["identity_snv"] == h["identity_snv"]]
        ident_snv = sum(i * w for i, w in snvw) / sum(w for _, w in snvw) if snvw else float("nan")
        frac = cov / mito_length
        # mtDNA content counts every nuclear base once: a nuclear copy of one tandem-repeat unit aligns to every copy of the unit
        # in the mitogenome, so mitogenome-side coverage alone can exceed the length of the locus itself
        nuc = sum(e - s for s, e in _union([(h["tstart"], h["tend"]) for h in hs]))
        content = min(cov, nuc)
        cls = "full-length" if content >= full_length * mito_length else ("large" if content >= large_bp else "fragment")
        rows.append({"chrom": chrom, "start": start + 1, "end": end, "span_bp": end - start, "n_hits": len(hs), "mito_bp": cov,
                     "nuclear_bp": nuc, "mtdna_bp": content,
                     "frac_of_mito": round(frac, 4), "identity": round(ident, 5), "identity_snv": round(ident_snv, 5),
                     "n_snv": sum(h["snv"] for h in hs), "n_indel": sum(h["indel"] for h in hs), "class": cls,
                     "strands": "".join(sorted(set(h["strand"] for h in hs))), "mito_intervals": ";".join(f"{s + 1}-{e}" for s, e in iv),
                     "intervals": iv})
    rows.sort(key=lambda r: (-r["mito_bp"], r["chrom"], r["start"]))
    for i, r in enumerate(rows, 1): r["locus"] = f"NUMT{i:03d}"
    return rows


def write_loci(rows: list[dict], path: Path) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(LOCUS_COLUMNS) + "\n")
        for r in rows: fh.write("\t".join(str(r[c]) for c in LOCUS_COLUMNS) + "\n")


def read_loci(path: Path) -> list[dict]:
    rows = []
    lines = Path(path).read_text().splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    for line in lines[1:]:
        if not line: continue
        r = dict(zip(header, line.split("\t")))
        for k in ("start", "end", "span_bp", "n_hits", "mito_bp", "n_snv", "n_indel"): r[k] = int(r[k])
        r["nuclear_bp"] = int(r["nuclear_bp"]) if r.get("nuclear_bp") not in (None, "") else r["span_bp"]
        r["mtdna_bp"] = int(r["mtdna_bp"]) if r.get("mtdna_bp") not in (None, "") else r["mito_bp"]
        for k in ("frac_of_mito", "identity", "identity_snv"): r[k] = float(r[k])
        r["intervals"] = [[int(a) - 1, int(b)] for a, b in (x.split("-") for x in r["mito_intervals"].split(";") if x)]
        rows.append(r)
    return rows


def mito_coverage(rows: list[dict], mito_length: int) -> list[tuple[int, int, float]]:
    """Per mitogenome column: number of NUMT loci covering it and the best identity among them."""
    n = [0] * mito_length; best = [0.0] * mito_length
    for r in rows:
        for s, e in r["intervals"]:
            for p in range(max(0, s), min(mito_length, e)):
                n[p] += 1; best[p] = max(best[p], r["identity"])
    return [(p + 1, n[p], round(best[p], 4)) for p in range(mito_length)]


def catalog_numts(mito_fasta: Path, nuclear_fasta: Path, outdir: Path, *, merge: int = 3000, use_blast: bool = False,
                  threads: int = 8) -> tuple[list[dict], dict]:
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    mito = read_fasta(mito_fasta); (mname, mseq), = mito.items(); L = len(mseq)
    hits = minimap2_hits(mito_fasta, nuclear_fasta, threads)
    if use_blast: hits += blastn_hits(mito_fasta, nuclear_fasta, outdir, threads)
    rows = merge_loci(hits, L, merge=merge)
    write_loci(rows, outdir / "numts.loci.tsv")
    cov = mito_coverage(rows, L)
    with open(outdir / "numts.mito_coverage.tsv", "w") as fh:
        fh.write("pos\tn_loci\tbest_identity\n")
        for p, n, b in cov: fh.write(f"{p}\t{n}\t{b}\n")
    summary = {
        "mito": mname, "mito_length": L, "n_loci": len(rows),
        "classes": dict(collections.Counter(r["class"] for r in rows)),
        "mito_bp_in_numts": sum(r["mito_bp"] for r in rows),
        "mtdna_bp_in_numts": sum(r["mtdna_bp"] for r in rows),
        "mito_columns_covered": sum(1 for _, n, _ in cov if n),
        "max_identity": max((r["identity"] for r in rows), default=None),
        "largest_locus_bp": max((r["mtdna_bp"] for r in rows), default=0),
        "hits": {"minimap2": sum(h["source"] == "minimap2" for h in hits), "blastn": sum(h["source"] == "blastn" for h in hits)},
    }
    (outdir / "numts.summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return rows, summary


# ------------------------------------------------------------------------------------------------ read classification
def _at_locus(loci_by_chrom, chrom, s, e, pad):
    """Name of a cataloged locus within ``pad`` bp of the half-open segment [s, e) (abutting intervals do not count)."""
    for ls, le, name in loci_by_chrom.get(chrom, []):
        if s < le + pad and e > ls - pad: return name
    return None


def _segment_vs_loci(loci_by_chrom, chrom, s, e, pad):
    """(locus name or None, bp of [s, e) inside that locus, bp outside it) for the nearest cataloged locus within ``pad``."""
    best = None
    for ls, le, name in loci_by_chrom.get(chrom, []):
        if s < le + pad and e > ls - pad:
            inside = max(0, min(e, le) - max(s, ls))
            if best is None or inside > best[1]:
                best = (name, inside, (e - s) - inside)
    return best or (None, 0, e - s)


def classify_reads(bam_path: Path, mito_name: str, loci: list[dict], *, pad: int = 2000, min_segment: int = 100,
                   min_flank: int = 500) -> list[dict]:
    """Alignment classes of reads mapped to nuclear + mitochondrial sequence together (primary + supplementary kept).

    These are placement classes, not proof of origin:
    ``mito_only`` / ``mito_multisegment`` (every segment on the mitogenome; several segments, e.g. a read crossing the origin);
    ``mito+nuclear_at_NUMT_locus`` (a nuclear segment near a cataloged locus with at least ``min_flank`` bp outside the locus:
    a candidate junction read with a real nuclear flank); ``mito+NUMT_homology`` (every nuclear segment lies mostly inside
    cataloged NUMT sequence: usually a mitochondrial read with a competing placement on the NUMT, no nuclear flank);
    ``mito+nuclear_elsewhere`` (other nuclear segments: an uncataloged or assembly-absent NUMT, or a library chimera);
    ``nuclear_only_at_NUMT_locus`` (every segment nuclear and overlapping a cataloged locus) and ``nuclear_only``.
    Nuclear-only reads can still carry mitochondrial sequence as unaligned insertions or clips (an assembly-absent NUMT)."""
    loci_by_chrom: dict[str, list] = collections.defaultdict(list)
    for r in loci: loci_by_chrom[r["chrom"]].append((r["start"] - 1, r["end"], r["locus"]))
    segs: dict[str, list] = collections.defaultdict(list); readlen: dict[str, int] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for a in bam.fetch(until_eof=True):
            if a.is_unmapped or a.is_secondary or a.query_alignment_length < min_segment: continue
            readlen[a.query_name] = max(readlen.get(a.query_name, 0), a.infer_read_length() or a.query_length)
            segs[a.query_name].append((a.reference_name, a.reference_start, a.reference_end, "-" if a.is_reverse else "+", a.query_alignment_length))
    out = []
    for name, ss in segs.items():
        mito = [s for s in ss if s[0] == mito_name]; nuc = [s for s in ss if s[0] != mito_name]
        if mito and not nuc:
            cls = "mito_only" if len(mito) == 1 else "mito_multisegment"; locus = ""
        elif nuc and not mito:
            labs = [_at_locus(loci_by_chrom, c, s, e, 0) for c, s, e, _, _ in nuc]          # must overlap, not just be near
            cls = "nuclear_only_at_NUMT_locus" if all(labs) else "nuclear_only"; locus = ",".join(sorted(set(x for x in labs if x)))
        else:
            calls = [_segment_vs_loci(loci_by_chrom, c, s, e, pad) for c, s, e, _, _ in nuc]
            flanks = [n for n, inside, outside in calls if n and outside >= min_flank]
            homology = [n for n, inside, outside in calls if n and inside >= 0.5 * (inside + outside)]
            if flanks:
                cls = "mito+nuclear_at_NUMT_locus"; locus = ",".join(sorted(set(flanks)))
            elif len(homology) == len(calls):
                cls = "mito+NUMT_homology"; locus = ",".join(sorted(set(homology)))
            else:
                cls = "mito+nuclear_elsewhere"; locus = ",".join(sorted(set(n for n, _, _ in calls if n)))
        detail = ";".join(f"{c}:{s}-{e}{st}" for c, s, e, st, _ in sorted(ss, key=lambda x: (x[0] != mito_name, x[1])))
        out.append({"read": name, "class": cls, "locus": locus, "read_length": readlen[name], "aligned_bp": sum(s[4] for s in ss), "segments": detail})
    return out


def write_read_classes(rows: list[dict], path: Path) -> dict:
    with open(path, "w") as fh:
        fh.write("read\tclass\tlocus\tread_length\taligned_bp\tsegments\n")
        for r in rows: fh.write("\t".join(str(r[k]) for k in ("read", "class", "locus", "read_length", "aligned_bp", "segments")) + "\n")
    return dict(collections.Counter(r["class"] for r in rows))


def tag_bam_with_classes(bam_in: Path, classes: list[dict], bam_out: Path) -> int:
    cls = {r["read"]: r["class"] for r in classes}
    inb = pysam.AlignmentFile(str(bam_in)); tmp = str(bam_out) + ".unsorted.bam"
    outb = pysam.AlignmentFile(tmp, "wb", template=inb); n = 0
    for a in inb:
        a.set_tag("PO", cls.get(a.query_name, "unclassified"), "Z"); outb.write(a); n += 1
    outb.close(); pysam.sort("-o", str(bam_out), tmp); pysam.index(str(bam_out)); Path(tmp).unlink()
    return n


def _aligned_bp(blocks, lo: int, hi: int) -> int:
    return sum(max(0, min(e, hi) - max(s, lo)) for s, e in blocks)


def junction_support(bam_path: Path, loci: list[dict], flank: int = 1000, min_frac: float = 0.9, min_mapq: int = 0) -> list[dict]:
    """Reads that cross a NUMT boundary with aligned bases on both sides: direct evidence that the insertion exists.

    For each boundary a read counts when at least ``min_frac`` of the ``flank`` bp outside the locus and of the first
    ``min(flank, locus length)`` bp inside it are covered by aligned bases of one alignment. Deletions and skips are not
    aligned bases, so a read from a haplotype without the insertion (a deletion spanning the locus) does not count.
    Secondary alignments are ignored; each read name counts once per boundary even if supplementary records repeat it."""
    out = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for r in loci:
            if r["chrom"] not in bam.references: continue
            chrom_len = bam.get_reference_length(r["chrom"])
            s0, e0 = r["start"] - 1, r["end"]
            inner = max(1, min(flank, e0 - s0))
            lo, hi = max(0, s0 - flank), min(chrom_len, e0 + flank)
            left, right = set(), set()
            for a in bam.fetch(r["chrom"], lo, hi):
                if a.is_unmapped or a.is_secondary or a.mapping_quality < min_mapq: continue
                blocks = a.get_blocks()
                if s0 > lo and _aligned_bp(blocks, lo, s0) >= min_frac * (s0 - lo) and _aligned_bp(blocks, s0, s0 + inner) >= min_frac * inner:
                    left.add(a.query_name)
                if hi > e0 and _aligned_bp(blocks, e0 - inner, e0) >= min_frac * inner and _aligned_bp(blocks, e0, hi) >= min_frac * (hi - e0):
                    right.add(a.query_name)
            out.append({"locus": r["locus"], "chrom": r["chrom"], "start": r["start"], "end": r["end"], "class": r["class"],
                        "left_junction_reads": len(left), "right_junction_reads": len(right)})
    return out


# ------------------------------------------------------------------------------------------------------------- figures
def _identity_color(ident: float, lo: float = 0.75):
    return plt.get_cmap("viridis")((max(lo, min(1.0, ident)) - lo) / (1 - lo))


def _identity_colorbar(fig, ax, *, fraction: float, label_size: float, tick_size: float):
    """Color bar for ``_identity_color``: identities below 75 % share the lowest color, so the bar is extended and the
    bottom tick reads "<=75"."""
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(75, 100))
    cb = fig.colorbar(sm, ax=ax, fraction=fraction, pad=0.01, extend="min")
    cb.set_ticks([75, 80, 85, 90, 95, 100]); cb.set_ticklabels(["≤75", "80", "85", "90", "95", "100"])
    cb.set_label("identity to mtDNA (%)", fontsize=label_size); cb.ax.tick_params(labelsize=tick_size)
    return cb


NUMT_CLASS_COLORS = {"full-length": "#d62728", "large": "#ff7f0e", "fragment": "#1f77b4"}


def class_labels(full_length: float = 0.95, large_bp: int = 5000) -> dict[str, str]:
    """Legend text that states what the locus classes mean: how much of the mitogenome one locus carries."""
    kb = f"{large_bp / 1000:g}"
    return {"full-length": f"full-length (\u2265 {full_length * 100:g} % of mtDNA)", "large": f"large (\u2265 {kb} kb of mtDNA)",
            "fragment": f"fragment (< {kb} kb of mtDNA)"}


def draw_numt_landscape(ax, loci: list[dict], chrom_lengths: dict[str, int], *, max_chroms: int = 40, min_chrom_bp: int = 0,
                        title: str | None = None, label_fontsize: float | None = None, min_visible_frac: float = 0.005) -> None:
    """Chromosomes as gray bars with every locus drawn at its true genomic span (box height = class, color = identity).
    A NUMT is a few kb on a chromosome of tens of Mb, so a locus narrower than ``min_visible_frac`` of the longest sequence
    is widened to that minimum, centered on the locus, and the title states the minimum and the true span range."""
    # show the long sequences (>= min_chrom_bp, default 2 % of the longest) plus every sequence that carries a locus
    if not min_chrom_bp and chrom_lengths:
        min_chrom_bp = int(0.02 * max(chrom_lengths.values()))
    with_loci = {r["chrom"] for r in loci}
    chroms = [c for c, l in sorted(chrom_lengths.items(), key=lambda kv: -kv[1]) if l >= min_chrom_bp or c in with_loci][:max_chroms]
    chroms = sorted(chroms, key=lambda c: -chrom_lengths[c])
    index = {c: i for i, c in enumerate(chroms)}
    if label_fontsize is None:
        label_fontsize = 6 if len(chroms) <= 28 else 5
    for c in chroms:
        ax.add_patch(Rectangle((0, index[c] - 0.3), chrom_lengths[c] / 1e6, 0.6, fc="#e5e5e5", ec="#999", lw=0.5))
    shown = 0
    longest = max((chrom_lengths[c] for c in chroms), default=1)
    min_w = min_visible_frac * longest
    widened = []
    for r in loci:
        if r["chrom"] not in index: continue
        shown += 1
        h = {"full-length": 0.62, "large": 0.46, "fragment": 0.3}[r["class"]]
        span = r["end"] - r["start"] + 1
        w = max(span, min_w)
        if span < min_w: widened.append(span)
        x0 = (r["start"] + r["end"]) / 2 - w / 2
        ax.add_patch(Rectangle((x0 / 1e6, index[r["chrom"]] - h / 2), w / 1e6, h, fc=_identity_color(r["identity"]),
                               ec="k" if r["class"] == "full-length" else "none", lw=0.8, zorder=3))
    ax.set_yticks(range(len(chroms))); ax.set_yticklabels(chroms, fontsize=label_fontsize)
    ax.set_ylim(-0.7, len(chroms) - 0.3); ax.set_xlim(0, max(chrom_lengths[c] for c in chroms) / 1e6 * 1.02 if chroms else 1)
    ax.set_xlabel("position (Mb)", fontsize=7); ax.tick_params(axis="x", labelsize=6); ax.invert_yaxis()
    n = collections.Counter(r["class"] for r in loci)
    note = "box height = class, color = identity to the mtDNA"
    if widened:
        lo, hi = min(widened) / 1000, max(widened) / 1000
        note += (f"; boxes \u2265 {min_w / 1000:.0f} kb wide (true spans {lo:.1f}-{hi:.1f} kb, "
                 f"up to {min_w / max(1, min(widened)):,.0f}x wider than real)")
    ax.set_title(title or f"NUMT landscape: {len(loci)} loci ({n.get('full-length', 0)} full-length, {n.get('large', 0)} large, "
                 f"{n.get('fragment', 0)} fragments); {shown} shown on {len(chroms)} sequences\n" + note, fontsize=6.5, loc="left")
    for s in ax.spines.values(): s.set_visible(False)


def draw_numt_sizes(ax, loci: list[dict], mito_length: int, *, full_length: float = 0.95, large_bp: int = 5000,
                    title: str | None = None, legend: bool = True, compact: bool = False) -> None:
    """Every locus as one bar of the mitochondrial sequence it carries, ranked, on an axis that spans the whole mitogenome.

    This is the panel that defines the classes: a bar reaching the dashed ``full-length`` line is a whole-mitogenome insertion, a
    bar past the ``large`` line carries >= ``large_bp`` of mtDNA (possibly as several pieces), everything shorter is a fragment.
    """
    rows = sorted(loci, key=lambda r: -r.get("mtdna_bp", r["mito_bp"]))
    n = len(rows); L = mito_length / 1000
    for i, r in enumerate(rows):
        ax.barh(i, r.get("mtdna_bp", r["mito_bp"]) / 1000, height=0.8 if n <= 60 else 1.0, color=NUMT_CLASS_COLORS[r["class"]], ec="none", zorder=3)
    labels = class_labels(full_length, large_bp)
    for x, txt, col in ((large_bp / 1000, "large", NUMT_CLASS_COLORS["large"]), (full_length * L, "full-length", NUMT_CLASS_COLORS["full-length"])):
        ax.axvline(x, color=col, lw=0.8, ls="--", zorder=2)
        ax.text(x, (n - 1) / 2, f"{txt} threshold ", color=col, fontsize=5.5, ha="right", va="center", rotation=90, zorder=4,
                bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.85))
    ax.set_xlim(0, L); ax.set_ylim(n - 0.4, -0.6)
    ax.set_yticks([]); ax.set_ylabel(f"{n} loci, ranked" if compact else f"{n} NUMT loci, ranked by mtDNA content", fontsize=6 if compact else 7)
    ax.set_xlabel("mtDNA in the locus (kb)", fontsize=7); ax.tick_params(axis="x", labelsize=6)
    top = ax.secondary_xaxis("top", functions=(lambda x: x / L * 100 if L else x, lambda p: p * L / 100 if L else p))
    top.set_xlabel("% of the mitogenome", fontsize=6); top.tick_params(labelsize=5)
    if legend:
        for cls, col in NUMT_CLASS_COLORS.items():
            ax.add_patch(Rectangle((0, 0), 0, 0, fc=col, ec="none", label=labels[cls]))
        ax.legend(fontsize=5 if compact else 5.5, loc="lower right", frameon=True, framealpha=0.9, edgecolor="none", handlelength=1.2)
    if title != "":
        ax.set_title(title or "mtDNA content per locus (class thresholds)", fontsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def add_annotation_strip(ax, gff: Path, mito_length: int, *, lane_in: float = 0.2, fontsize: float = 5.0):
    """The gene annotation as its own labeled track directly under an axes in mitogenome coordinates (kb).

    The strip takes the bottom of ``ax``'s box and is attached to ``ax`` (it moves and resizes with it) with a shared x axis.
    Lanes: protein-coding and rRNA genes, named where the name fits inside the block, then tRNA genes, each split into + and -
    strand lanes when both strands occur. Features that run past the end of the molecule are folded onto its start.
    Unannotated stretches of at least 2 % of the molecule are marked "non-coding". The x tick labels and axis label move to the
    strip. Returns the strip axes, or ``None`` when the GFF has no features. Colors are redwood's annotation colors."""
    from .renderer import FEATURE_COLORS, parse_gff
    feats = parse_gff(Path(gff))
    if not feats:
        return None
    pieces = []                                              # (start, stop, feature): 0-based half-open, folded at the origin
    for f in feats:
        s0, e0 = int(f["start"]), int(f["stop"])
        if e0 > mito_length:
            pieces.append((s0, mito_length, f)); pieces.append((0, e0 - mito_length, f))
        else:
            pieces.append((s0, e0, f))
    kinds = []
    for kind, is_trna in (("genes", False), ("tRNAs", True)):
        strands = sorted({str(f.get("strand", "+")) for _, _, f in pieces if (str(f["type"]) == "tRNA") == is_trna})
        for st in strands:
            kinds.append((kind if len(strands) == 1 else f"{kind} ({st})", is_trna, st))
    n = len(kinds)
    fig = ax.figure
    box = ax.get_position(original=True)
    frac = min(0.4 * box.height, n * lane_in / fig.get_size_inches()[1])
    ax.set_position([box.x0, box.y0 + frac, box.width, box.height - frac])
    rel = frac / (box.height - frac)
    ann = ax.inset_axes([0, -rel, 1, rel], sharex=ax)
    ax.tick_params(labelbottom=False); xlabel = ax.get_xlabel(); ax.set_xlabel("")
    kb_per_pt = (mito_length / 1000) / (box.width * fig.get_size_inches()[0] * 72)
    lane_y = {(is_trna, st): n - 1 - i for i, (_, is_trna, st) in enumerate(kinds)}   # first lane on top
    for s0, e0, f in pieces:
        is_trna = str(f["type"]) == "tRNA"
        y = lane_y[(is_trna, str(f.get("strand", "+")))]
        start_kb, stop_kb = s0 / 1000, e0 / 1000
        color = FEATURE_COLORS.get(str(f["type"]), "#d08c35")
        h = 0.7 if is_trna else 0.8
        ann.add_patch(Rectangle((start_kb, y + (1 - h) / 2), max(stop_kb - start_kb, 0.5 * kb_per_pt), h, fc=color, ec="none"))
        name = str(f["name"])
        if not is_trna and (stop_kb - start_kb) >= 1.15 * len(name) * fontsize * 0.6 * kb_per_pt:
            ann.text((start_kb + stop_kb) / 2, y + 0.5, name, ha="center", va="center", fontsize=fontsize, color="white",
                     fontweight="bold", clip_on=True)
    first_gene_lane = next((lane_y[(t, st)] for _, t, st in kinds if not t), n - 1)
    end = 0
    for s0, e0 in sorted((s0, e0) for s0, e0, _ in pieces) + [(mito_length, mito_length)]:
        if s0 - end >= 0.02 * mito_length and (s0 - end) / 1000 >= 1.15 * 10 * fontsize * 0.6 * kb_per_pt:
            ann.text((s0 + end) / 2000, first_gene_lane + 0.5, "non-coding", ha="center", va="center", fontsize=fontsize,
                     color="#667085", fontstyle="italic")
        end = max(end, e0)
    ticks = sorted((lane_y[(t, st)] + 0.5, k) for k, t, st in kinds)
    ann.set_ylim(0, n); ann.set_yticks([y for y, _ in ticks]); ann.set_yticklabels([k for _, k in ticks], fontsize=5)
    ann.set_ylabel("annotation", fontsize=5.5)
    ann.tick_params(axis="y", length=0); ann.tick_params(axis="x", labelsize=6)
    ann.set_xlabel(xlabel, fontsize=7)
    for side in ("top", "right"):
        ann.spines[side].set_visible(False)
    return ann


def draw_numt_catalog(ax, loci: list[dict], mito_length: int, gff: Path | None = None, *, title: str | None = None,
                      full_length: float = 0.95, large_bp: int = 5000, legend_ncol: int = 3) -> None:
    """Each locus as horizontal segments over the mitogenome coordinates it covers, at y = identity, colored by class.

    A locus whose alignments cover the mitogenome in several intervals keeps one identity and one class; the intervals are
    joined by a dotted line (in the legend as "alignments in one locus"). They can be separate pieces of the insertion or
    alternative placements of the same nuclear bases on repeated mitochondrial sequence. The class is the mitochondrial
    content of the whole locus, counted once per nuclear base. With ``gff`` the annotation is drawn as a labeled strip under the axes (:func:`add_annotation_strip`).
    """
    from matplotlib.lines import Line2D
    cc = NUMT_CLASS_COLORS
    for r in loci:
        y = r["identity"] * 100; iv = sorted(r["intervals"])
        if len(iv) > 1:
            ax.plot([iv[0][0] / 1000, iv[-1][1] / 1000], [y, y], color=cc[r["class"]], lw=0.6, ls=":", alpha=0.7, zorder=2)
        for s, e in iv:
            ax.plot([s / 1000, e / 1000], [y, y], color=cc[r["class"]], lw=2.4 if r["class"] != "fragment" else 1.2,
                    solid_capstyle="butt", alpha=0.9, zorder=3)
    ax.set_xlim(0, mito_length / 1000); ax.set_ylabel("identity to mtDNA (%)", fontsize=7); ax.set_xlabel("mitogenome position (kb)", fontsize=7)
    ax.tick_params(labelsize=6)
    lo = min([r["identity"] * 100 for r in loci] + [95]); bottom = max(60, lo - 3)
    ax.set_ylim(bottom - 0.24 * (100.5 - bottom), 100.5)          # an empty band under the data holds the legend
    ax.set_yticks([t for t in ax.get_yticks() if bottom - 1e-9 <= t <= 100])
    kb = f"{large_bp / 1000:g}"
    labels = {"full-length": f"full-length (\u2265{full_length * 100:g}%)", "large": f"large (\u2265{kb} kb)",
              "fragment": f"fragment (<{kb} kb)"}                     # short: panel b's legend spells the classes out
    handles = [Line2D([], [], color=col, lw=2, label=labels[cls]) for cls, col in cc.items()]
    handles.append(Line2D([], [], color="#555555", lw=0.9, ls=":", label="alignments in one locus"))
    ax.legend(handles=handles, fontsize=6 if legend_ncol > 1 else 5.5, loc="lower right", frameon=True, framealpha=0.9, edgecolor="none",
              ncol=4 if legend_ncol > 1 else 1, columnspacing=1.0, handlelength=1.8)
    ax.set_title(title or f"NUMT catalog: {len(loci)} loci, {sum(r.get('mtdna_bp', r['mito_bp']) for r in loci):,} bp of mtDNA in the nuclear genome\n"
                 "class = mtDNA content of the whole locus", fontsize=6.5 if legend_ncol > 1 else 6, loc="left")
    if gff is not None and Path(gff).exists():
        add_annotation_strip(ax, Path(gff), mito_length)


def plot_numt_figures(loci: list[dict], mito_length: int, chrom_lengths: dict[str, int], outdir: Path, gff: Path | None = None,
                      dpi: int = 200) -> list[Path]:
    outdir = Path(outdir); outs = []
    fig, ax = plt.subplots(figsize=(7.5, max(2.5, 0.16 * min(40, len(chrom_lengths)) + 1)))
    draw_numt_landscape(ax, loci, chrom_lengths)
    _identity_colorbar(fig, ax, fraction=0.03, label_size=7, tick_size=6)
    fig.savefig(outdir / "numts.landscape.png", dpi=dpi, bbox_inches="tight"); fig.savefig(outdir / "numts.landscape.pdf", bbox_inches="tight"); plt.close(fig)
    outs.append(outdir / "numts.landscape.png")
    fig, ax = plt.subplots(figsize=(7.5, 3)); draw_numt_catalog(ax, loci, mito_length, gff)
    fig.savefig(outdir / "numts.catalog.png", dpi=dpi, bbox_inches="tight"); fig.savefig(outdir / "numts.catalog.pdf", bbox_inches="tight"); plt.close(fig)
    outs.append(outdir / "numts.catalog.png")
    fig, ax = plt.subplots(figsize=(3.2, max(2.5, min(6.0, 0.05 * len(loci) + 1.2)))); draw_numt_sizes(ax, loci, mito_length)
    fig.savefig(outdir / "numts.sizes.png", dpi=dpi, bbox_inches="tight"); fig.savefig(outdir / "numts.sizes.pdf", bbox_inches="tight"); plt.close(fig)
    outs.append(outdir / "numts.sizes.png")
    return outs


def plot_composite_figure(*, mito_fasta: Path, loci: list[dict], chrom_lengths: dict[str, int], out_base: Path, gff: Path | None = None,
                          main_bam: Path | None = None, rnaseq_bam: Path | None = None, variant_table: Path | None = None,
                          title: str | None = None, dpi: int = 300, fileforms=("png", "pdf"), page: tuple[float, float] = (6.5, 9.0),
                          **circular_kwargs) -> list[Path]:
    """One page per mitogenome: the redwood circular map with its ring key (a) and, beside it, the per-locus mtDNA-content panel
    that defines the NUMT classes (b); then the NUMT landscape (c) and the NUMT catalog (d) across the full width.

    The page is exactly ``page`` inches (default 6.5 x 9, the text area of letter paper with 1-inch margins, so the file prints on
    letter or A4 and drops into a Word / Google document at full text width without scaling); margins are laid out inside the
    figure and the output is not cropped, so the PDF page and the PNG pixel size are fixed."""
    from .renderer import draw_circular_plot, read_reference
    reference = read_reference(Path(mito_fasta)); L = len(reference)
    w, h = page
    left, right, top, bottom = 0.105, 0.925, 0.965, 0.06
    nchrom = len({r["chrom"] for r in loci} | {c for c, l in chrom_lengths.items() if l >= 0.02 * max(chrom_lengths.values())}) if chrom_lengths else 1
    nchrom = min(40, nchrom)
    h_map = 0.37 * h
    h_land = min(0.26 * h, max(0.15 * h, (0.085 * nchrom + 0.5) * h / 9.0))
    h_cat = 0.23 * h
    fig = plt.figure(figsize=(w, h))
    gs = fig.add_gridspec(3, 2, width_ratios=[0.64, 0.36], height_ratios=[h_map, h_land, h_cat], hspace=0.40, wspace=0.30,
                          left=left, right=right, top=top, bottom=bottom)
    ax0 = fig.add_subplot(gs[0, 0]); ax_sz = fig.add_subplot(gs[0, 1]); ax1 = fig.add_subplot(gs[1, :]); ax2 = fig.add_subplot(gs[2, :])
    ax0.set_anchor("N")                                     # map and class panel share their top edge
    draw_circular_plot(ax0, length=L, reference=reference, gff=Path(gff) if gff else None, main_bam=Path(main_bam) if main_bam else None,
                       rnaseq_bam=Path(rnaseq_bam) if rnaseq_bam else None, variant_table=Path(variant_table) if variant_table else None,
                       title=title, **circular_kwargs)
    for s in ax0.spines.values(): s.set_visible(False)
    draw_numt_sizes(ax_sz, loci, L, title="")
    row_pt = ax1.get_position().height * h * 72 / (nchrom + 0.4)      # points available per chromosome row
    draw_numt_landscape(ax1, loci, chrom_lengths, label_fontsize=max(5.0, min(6.0, 0.8 * row_pt)))
    _identity_colorbar(fig, ax1, fraction=0.02, label_size=6, tick_size=5)
    draw_numt_catalog(ax2, loci, L, gff)
    for letter, ax in zip("abcd", (ax0, ax_sz, ax1, ax2)):        # panel letters in the page margin, clear of titles and tick labels
        box = ax.get_position()
        x = max(0.012, box.x0 - 0.085) if ax is not ax_sz else box.x0 - 0.085
        fig.text(x, min(0.985, box.y1 + 0.006), letter, fontsize=10, fontweight="bold", va="bottom", ha="left")
    outs = []
    for ff in fileforms:
        p = Path(f"{out_base}.{ff}"); fig.savefig(p, dpi=dpi); outs.append(p)
    plt.close(fig)
    return outs


# --------------------------------------------------------------------------------------------------------------- CLI
def run_numts(args) -> dict:
    from .workflow import require_tool
    require_tool("minimap2")
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    loci, summary = catalog_numts(Path(args.mito_fasta), Path(args.nuclear_fasta), outdir, merge=args.merge, use_blast=args.blastn, threads=args.threads)
    chrom_lengths = {}
    fai = Path(str(args.nuclear_fasta) + ".fai")
    if not fai.exists(): pysam.faidx(str(args.nuclear_fasta))
    for line in fai.read_text().splitlines():
        f = line.split("\t"); chrom_lengths[f[0]] = int(f[1])
    mito_name = summary["mito"]
    # reads: map to nuclear + mito together (-Y, supplementary kept) and classify
    bam = Path(args.bam) if args.bam else None
    if args.long_reads:
        require_tool("samtools")
        combined = outdir / "nuclear_plus_mito.fa"
        with open(combined, "w") as fh:
            for src in (args.nuclear_fasta, args.mito_fasta):
                fh.write(Path(src).read_text().rstrip("\n") + "\n")
        bam = outdir / "reads_vs_nuclear_plus_mito.bam"
        cmd = ["minimap2", "-t", str(args.threads), "-a", "-x", args.long_read_preset, "-Y", "--secondary=no", str(combined)] + [str(p) for p in args.long_reads]
        with open(outdir / "map.log", "w") as log:
            p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log)
            p2 = subprocess.Popen(["samtools", "view", "-b", "-F", "0x104", "-"], stdin=p1.stdout, stdout=subprocess.PIPE)
            p1.stdout.close()                     # so minimap2 sees a broken pipe if samtools exits early
            p3 = subprocess.run(["samtools", "sort", "-@", "4", "-o", str(bam), "-"], stdin=p2.stdout)
            p2.stdout.close()
            status = {"minimap2": p1.wait(), "samtools view": p2.wait(), "samtools sort": p3.returncode}
            failed = {k: v for k, v in status.items() if v != 0}
            if failed:
                raise RuntimeError(f"read mapping failed ({failed}); see {outdir / 'map.log'}")
        subprocess.run(["samtools", "index", str(bam)], check=True)
    if bam is not None:
        classes = classify_reads(bam, mito_name, loci, pad=args.pad, min_flank=getattr(args, "min_flank", 500))
        summary["read_classes"] = write_read_classes(classes, outdir / "numts.read_classes.tsv")
        js = junction_support(bam, loci, flank=args.flank)
        with open(outdir / "numts.junction_support.tsv", "w") as fh:
            fh.write("locus\tchrom\tstart\tend\tclass\tleft_junction_reads\tright_junction_reads\n")
            for r in js: fh.write("\t".join(str(r[k]) for k in ("locus", "chrom", "start", "end", "class", "left_junction_reads", "right_junction_reads")) + "\n")
        summary["loci_with_both_junctions_supported"] = sum(1 for r in js if r["left_junction_reads"] and r["right_junction_reads"])
        if args.mito_bam:
            summary["tagged_records"] = tag_bam_with_classes(Path(args.mito_bam), classes, outdir / "reads_vs_mito.PO.bam")
    figs = plot_numt_figures(loci, summary["mito_length"], chrom_lengths, outdir, gff=Path(args.gff) if args.gff else None, dpi=args.dpi)
    summary["figures"] = [str(p) for p in figs]
    if args.figure:
        outs = plot_composite_figure(mito_fasta=Path(args.mito_fasta), loci=loci, chrom_lengths=chrom_lengths, out_base=outdir / args.figure_name,
                                     gff=Path(args.gff) if args.gff else None, main_bam=Path(args.figure_bam) if args.figure_bam else None,
                                     rnaseq_bam=Path(args.figure_rnaseq_bam) if args.figure_rnaseq_bam else None,
                                     variant_table=Path(args.figure_variant_table) if args.figure_variant_table else None, dpi=args.dpi)
        summary["composite_figure"] = [str(p) for p in outs]
    (outdir / "numts.summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary
