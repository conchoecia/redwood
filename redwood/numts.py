"""NUMT catalogue, read classification and figures.

A nuclear genome usually carries copies of mitochondrial sequence (NUMTs). They confuse mitogenome assembly (an assembler may
extend the mitochondrion into a NUMT flank, or build the NUMT instead of the molecule) and read-population analyses (NUMT-derived
reads look like a divergent mitochondrial haplotype). This module

* catalogues NUMT loci: the mitogenome is aligned to the nuclear assembly (minimap2 ``asm20``; optionally ``blastn -task
  dc-megablast`` for short diverged fragments), hits on one chromosome within ``merge`` bp are merged into a locus, and each locus
  gets its span, the mitogenome intervals it covers, the fraction of the molecule covered, SNV-only and total identity, and a class
  (``full-length`` >= 95 % of the molecule, ``large`` >= 5 kb, ``fragment``);
* classifies long reads mapped to nuclear + mitochondrial sequence together (``minimap2 -Y --secondary=no``, supplementary
  alignments kept) by the origin of their segments: ``mito_only``, ``mito_multisegment`` (two mito segments, e.g. a read crossing
  the origin of a circular molecule), ``mito+nuclear_at_NUMT_locus`` (a junction read of a catalogued NUMT), ``mito+nuclear_elsewhere``
  (an uncatalogued NUMT or a library chimera), ``nuclear_only_at_NUMT_locus`` (a NUMT read whose whole length is explained by the
  nuclear locus) and ``nuclear_only``;
* counts reads spanning each nuclear-mitochondrial junction in a single alignment (junction support);
* draws the NUMT landscape (chromosomes with loci coloured by identity), the NUMT-versus-mitogenome catalogue (which part of the
  molecule each locus covers, at what identity) and a composite figure (circular map on top, landscape and catalogue below).
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
}
LOCUS_COLUMNS = ["locus", "chrom", "start", "end", "span_bp", "n_hits", "mito_bp", "frac_of_mito", "identity", "identity_snv",
                 "n_snv", "n_indel", "class", "strands", "mito_intervals"]


# ----------------------------------------------------------------------------------------------------------- catalogue
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


def blastn_hits(mito_fasta: Path, nuclear_fasta: Path, workdir: Path, threads: int = 8, min_len: int = 150) -> list[dict]:
    """Short / diverged fragments that asm20 chaining misses (dc-megablast). Needs BLAST+ on PATH."""
    if not (shutil.which("blastn") and shutil.which("makeblastdb")):
        return []
    db = workdir / "nuclear_db"
    subprocess.run(["makeblastdb", "-in", str(nuclear_fasta), "-dbtype", "nucl", "-out", str(db)], check=True, capture_output=True)
    out = subprocess.run(["blastn", "-task", "dc-megablast", "-query", str(mito_fasta), "-db", str(db), "-evalue", "1e-10",
                          "-num_threads", str(threads), "-max_target_seqs", "50", "-max_hsps", "5000",
                          "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"],
                         capture_output=True, text=True, check=True).stdout
    hits = []
    for line in out.splitlines():
        q, t, pid, ln, mm, go, qs, qe, ss, se, ev, bs = line.split("\t")
        if int(ln) < min_len: continue
        ss, se = int(ss), int(se); strand = "+" if se >= ss else "-"
        hits.append({"chrom": t, "tstart": min(ss, se) - 1, "tend": max(ss, se), "qstart": int(qs) - 1, "qend": int(qe), "strand": strand,
                     "alen": int(ln), "identity": float(pid) / 100, "identity_snv": (int(ln) - int(mm)) / int(ln), "snv": int(mm),
                     "indel": int(go), "source": "blastn"})
    return hits


def _union(intervals: list[tuple[int, int]]) -> list[list[int]]:
    out: list[list[int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out


def merge_loci(hits: list[dict], mito_length: int, merge: int = 3000, full_length: float = 0.95, large_bp: int = 5000) -> list[dict]:
    """Cluster hits per chromosome (any hit within ``merge`` bp of the running locus joins it) and summarise each locus."""
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
        cls = "full-length" if frac >= full_length else ("large" if cov >= large_bp else "fragment")
        rows.append({"chrom": chrom, "start": start + 1, "end": end, "span_bp": end - start, "n_hits": len(hs), "mito_bp": cov,
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
    header = lines[0].split("\t")
    for line in lines[1:]:
        if not line: continue
        r = dict(zip(header, line.split("\t")))
        for k in ("start", "end", "span_bp", "n_hits", "mito_bp", "n_snv", "n_indel"): r[k] = int(r[k])
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
        "mito_columns_covered": sum(1 for _, n, _ in cov if n),
        "max_identity": max((r["identity"] for r in rows), default=None),
        "largest_locus_bp": max((r["mito_bp"] for r in rows), default=0),
        "hits": {"minimap2": sum(h["source"] == "minimap2" for h in hits), "blastn": sum(h["source"] == "blastn" for h in hits)},
    }
    (outdir / "numts.summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return rows, summary


# ------------------------------------------------------------------------------------------------ read classification
def _at_locus(loci_by_chrom, chrom, s, e, pad):
    for ls, le, name in loci_by_chrom.get(chrom, []):
        if s <= le + pad and e >= ls - pad: return name
    return None


def classify_reads(bam_path: Path, mito_name: str, loci: list[dict], *, pad: int = 2000, min_segment: int = 100) -> list[dict]:
    """Segment classes of reads mapped to nuclear + mitochondrial sequence together (primary + supplementary kept)."""
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
        if mito and not nuc: cls = "mito_only" if len(mito) == 1 else "mito_multisegment"; locus = ""
        elif nuc and not mito:
            labs = [_at_locus(loci_by_chrom, c, s, e, pad) for c, s, e, _, _ in nuc]
            cls = "nuclear_only_at_NUMT_locus" if all(labs) else "nuclear_only"; locus = ",".join(sorted(set(x for x in labs if x)))
        else:
            labs = [_at_locus(loci_by_chrom, c, s, e, pad) for c, s, e, _, _ in nuc]
            cls = "mito+nuclear_at_NUMT_locus" if all(labs) else "mito+nuclear_elsewhere"; locus = ",".join(sorted(set(x for x in labs if x)))
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


def junction_support(bam_path: Path, loci: list[dict], flank: int = 1000) -> list[dict]:
    """Reads whose single alignment spans a NUMT junction (locus boundary +- flank): direct evidence that the insertion exists in the animal."""
    out = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for r in loci:
            if r["chrom"] not in bam.references: continue
            left = right = 0
            for a in bam.fetch(r["chrom"], max(0, r["start"] - flank - 1), r["end"] + flank):
                if a.is_unmapped or a.is_secondary: continue
                if a.reference_start <= r["start"] - 1 - flank and a.reference_end >= r["start"] - 1 + flank: left += 1
                if a.reference_start <= r["end"] - flank and a.reference_end >= r["end"] + flank: right += 1
            out.append({"locus": r["locus"], "chrom": r["chrom"], "start": r["start"], "end": r["end"], "class": r["class"],
                        "left_junction_reads": left, "right_junction_reads": right})
    return out


# ------------------------------------------------------------------------------------------------------------- figures
def _identity_color(ident: float, lo: float = 0.75):
    return plt.get_cmap("viridis")((max(lo, min(1.0, ident)) - lo) / (1 - lo))


def draw_numt_landscape(ax, loci: list[dict], chrom_lengths: dict[str, int], *, max_chroms: int = 40, min_chrom_bp: int = 0,
                        title: str | None = None) -> None:
    # show the long sequences (>= min_chrom_bp, default 2 % of the longest) plus every sequence that carries a locus
    if not min_chrom_bp and chrom_lengths:
        min_chrom_bp = int(0.02 * max(chrom_lengths.values()))
    with_loci = {r["chrom"] for r in loci}
    chroms = [c for c, l in sorted(chrom_lengths.items(), key=lambda kv: -kv[1]) if l >= min_chrom_bp or c in with_loci][:max_chroms]
    chroms = sorted(chroms, key=lambda c: -chrom_lengths[c])
    index = {c: i for i, c in enumerate(chroms)}
    for c in chroms:
        ax.add_patch(Rectangle((0, index[c] - 0.3), chrom_lengths[c] / 1e6, 0.6, fc="#e5e5e5", ec="#999", lw=0.5))
    shown = 0
    for r in loci:
        if r["chrom"] not in index: continue
        shown += 1
        h = {"full-length": 0.62, "large": 0.46, "fragment": 0.3}[r["class"]]
        w = max(0.3, r["span_bp"] / 1e6)
        ax.add_patch(Rectangle((r["start"] / 1e6 - 0.15, index[r["chrom"]] - h / 2), w + 0.15, h, fc=_identity_color(r["identity"]),
                               ec="k" if r["class"] == "full-length" else "none", lw=0.8, zorder=3))
    ax.set_yticks(range(len(chroms))); ax.set_yticklabels(chroms, fontsize=6)
    ax.set_ylim(-0.7, len(chroms) - 0.3); ax.set_xlim(0, max(chrom_lengths[c] for c in chroms) / 1e6 * 1.02 if chroms else 1)
    ax.set_xlabel("position (Mb)", fontsize=7); ax.tick_params(axis="x", labelsize=6); ax.invert_yaxis()
    n = collections.Counter(r["class"] for r in loci)
    ax.set_title(title or f"NUMT landscape: {len(loci)} loci ({n.get('full-length', 0)} full-length, {n.get('large', 0)} >= 5 kb, "
                 f"{n.get('fragment', 0)} fragments); {shown} shown on {len(chroms)} sequences", fontsize=7)
    for s in ax.spines.values(): s.set_visible(False)


def draw_numt_catalog(ax, loci: list[dict], mito_length: int, gff: Path | None = None, *, title: str | None = None) -> None:
    """Each locus as horizontal segments over the mitogenome coordinates it covers, at y = identity, coloured by class."""
    cc = {"full-length": "#d62728", "large": "#ff7f0e", "fragment": "#1f77b4"}
    for r in loci:
        for s, e in r["intervals"]:
            ax.plot([s / 1000, e / 1000], [r["identity"] * 100] * 2, color=cc[r["class"]], lw=2.2 if r["class"] != "fragment" else 1.2,
                    solid_capstyle="butt", alpha=0.9)
    ax.set_xlim(0, mito_length / 1000); ax.set_ylabel("identity to mtDNA (%)", fontsize=7); ax.set_xlabel("mitogenome position (kb)", fontsize=7)
    ax.tick_params(labelsize=6)
    lo = min([r["identity"] * 100 for r in loci] + [95]); ax.set_ylim(max(60, lo - 3), 100.5)
    if gff is not None and Path(gff).exists():
        from .renderer import parse_gff, FEATURE_COLORS
        y0 = ax.get_ylim()[0]
        for f in parse_gff(Path(gff)):
            ax.add_patch(Rectangle((int(f["start"]) / 1000, y0), (int(f["stop"]) - int(f["start"])) / 1000, 1.2,
                                   fc=FEATURE_COLORS.get(str(f["type"]), "#d08c35"), ec="none", clip_on=False))
    for cls, col in cc.items():
        ax.plot([], [], color=col, lw=2, label=cls)
    ax.legend(fontsize=6, loc="lower right", frameon=False, ncol=3)
    n = collections.Counter(r["class"] for r in loci)
    ax.set_title(title or f"NUMT catalogue: {len(loci)} loci, {sum(r['mito_bp'] for r in loci):,} bp of mitochondrial sequence in the nuclear genome", fontsize=7)


def plot_numt_figures(loci: list[dict], mito_length: int, chrom_lengths: dict[str, int], outdir: Path, gff: Path | None = None,
                      dpi: int = 200) -> list[Path]:
    outdir = Path(outdir); outs = []
    fig, ax = plt.subplots(figsize=(7.5, max(2.5, 0.16 * min(40, len(chrom_lengths)) + 1)))
    draw_numt_landscape(ax, loci, chrom_lengths)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(75, 100)); cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.01)
    cb.set_label("identity to mtDNA (%)", fontsize=7); cb.ax.tick_params(labelsize=6)
    fig.savefig(outdir / "numts.landscape.png", dpi=dpi, bbox_inches="tight"); fig.savefig(outdir / "numts.landscape.pdf", bbox_inches="tight"); plt.close(fig)
    outs.append(outdir / "numts.landscape.png")
    fig, ax = plt.subplots(figsize=(7.5, 3)); draw_numt_catalog(ax, loci, mito_length, gff)
    fig.savefig(outdir / "numts.catalog.png", dpi=dpi, bbox_inches="tight"); fig.savefig(outdir / "numts.catalog.pdf", bbox_inches="tight"); plt.close(fig)
    outs.append(outdir / "numts.catalog.png")
    return outs


def plot_composite_figure(*, mito_fasta: Path, loci: list[dict], chrom_lengths: dict[str, int], out_base: Path, gff: Path | None = None,
                          main_bam: Path | None = None, rnaseq_bam: Path | None = None, variant_table: Path | None = None,
                          title: str | None = None, dpi: int = 300, fileforms=("png", "pdf"), **circular_kwargs) -> list[Path]:
    """One figure per mitogenome: redwood circular map on top, NUMT landscape and NUMT catalogue below."""
    from .renderer import draw_circular_plot, read_reference
    reference = read_reference(Path(mito_fasta)); L = len(reference)
    nchrom = min(40, len(chrom_lengths)); h_land = max(2.2, 0.14 * nchrom + 0.8)
    fig = plt.figure(figsize=(7.2, 7.2 + h_land + 2.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[7.2, h_land, 2.4], hspace=0.25)
    ax0 = fig.add_subplot(gs[0]); ax1 = fig.add_subplot(gs[1]); ax2 = fig.add_subplot(gs[2])
    draw_circular_plot(ax0, length=L, reference=reference, gff=Path(gff) if gff else None, main_bam=Path(main_bam) if main_bam else None,
                       rnaseq_bam=Path(rnaseq_bam) if rnaseq_bam else None, variant_table=Path(variant_table) if variant_table else None,
                       title=title, **circular_kwargs)
    for s in ax0.spines.values(): s.set_visible(False)
    draw_numt_landscape(ax1, loci, chrom_lengths)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(75, 100)); cb = fig.colorbar(sm, ax=ax1, fraction=0.03, pad=0.01)
    cb.set_label("identity to mtDNA (%)", fontsize=6); cb.ax.tick_params(labelsize=5)
    draw_numt_catalog(ax2, loci, L, gff)
    for i, ax in enumerate((ax0, ax1, ax2)):
        ax.text(-0.02, 1.0, "abc"[i], transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")
    outs = []
    for ff in fileforms:
        p = Path(f"{out_base}.{ff}"); fig.savefig(p, dpi=dpi, bbox_inches="tight"); outs.append(p)
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
            subprocess.run(["samtools", "sort", "-@", "4", "-o", str(bam), "-"], stdin=p2.stdout, check=True)
        subprocess.run(["samtools", "index", str(bam)], check=True)
    if bam is not None:
        classes = classify_reads(bam, mito_name, loci, pad=args.pad)
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
