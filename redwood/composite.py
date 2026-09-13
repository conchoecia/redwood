"""Journal-style composite figures.

``redwood composite`` draws the circular mitogenome map together with the NUMT panels (mtDNA content per locus, landscape,
catalog) as one figure laid out to a journal's geometry: figure width, figure text sizes, panel labels and, on a full page, the
legend under the figure. Styles are presets of measured values that every option can override.

The ``nature-communications`` style was measured from 2025 Nature Communications article PDFs: 210 x 279 mm pages; figures
placed about 170 mm wide, centered on the 183.8 mm two-column measure (two 90 mm columns, 3.9 mm gutter), starting 17.2 mm
below the page top; legends 3 mm below the figure in 7 pt type on 10 pt leading, in two columns; figure text 5-7 pt sans
serif; panel labels 8 pt bold lowercase. A 170 mm figure also matches the 6.7 in text width of a letter-paper supplement.

Legend text uses light markup: ``**bold**``, ``*italic*``, ``***bold italic***``; ``{field}`` placeholders are filled from the
data ({length}, {n_drawn}, {n_reads}, {min_box_kb}, {span_lo_kb}, {span_hi_kb}, {max_factor}), the panel letters ({a} .. {d},
bold) and ``--field KEY=VALUE``.
"""
from __future__ import annotations

import glob
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

MM = 1 / 25.4
PAGE_SIZES_MM = {"letter": (215.9, 279.4), "a4": (210.0, 297.0), "nature-communications": (210.0, 279.0)}
SANS_SERIF = ("Helvetica", "Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans")
DEFAULT_TITLE = "Mitochondrial genome and nuclear mitochondrial insertions (NUMTs)"


@dataclass(frozen=True)
class JournalStyle:
    """Page geometry and typography of one journal (lengths in mm, sizes in pt)."""
    name: str
    page_mm: tuple[float, float]
    figure_width_mm: float
    measure_mm: float
    gutter_mm: float
    top_mm: float
    bottom_mm: float              # text block bottom, measured from the page bottom
    legend_gap_mm: float
    legend_pt: float
    leading_pt: float
    legend_columns: int
    text_pt: tuple[float, float]
    panel_label_pt: float
    panel_label_case: str         # "lower" or "upper"
    fonts: tuple[str, ...]
    figure_prefix: str
    supplementary_prefix: str
    label_separator: str


STYLES = {
    "nature-communications": JournalStyle(
        name="nature-communications", page_mm=(210.0, 279.0), figure_width_mm=170.0, measure_mm=183.8, gutter_mm=3.9,
        top_mm=17.2, bottom_mm=19.9, legend_gap_mm=3.0, legend_pt=7.0, leading_pt=10.0, legend_columns=2, text_pt=(5.0, 7.0),
        panel_label_pt=8.0, panel_label_case="lower", fonts=SANS_SERIF, figure_prefix="Fig.",
        supplementary_prefix="Supplementary Fig.", label_separator=" | "),
    "generic": JournalStyle(
        name="generic", page_mm=(215.9, 279.4), figure_width_mm=165.1, measure_mm=165.1, gutter_mm=5.0, top_mm=25.4,
        bottom_mm=25.4, legend_gap_mm=4.0, legend_pt=8.0, leading_pt=10.0, legend_columns=1, text_pt=(5.0, 8.0),
        panel_label_pt=10.0, panel_label_case="upper", fonts=SANS_SERIF, figure_prefix="Figure",
        supplementary_prefix="Supplementary Figure", label_separator=". "),
}


# ------------------------------------------------------------------------------------------------------------ helpers
def parse_length_mm(value) -> float:
    """'170', '170mm', '17cm' or '6.7in' -> millimeters."""
    if isinstance(value, (int, float)):
        return float(value)
    m = re.fullmatch(r"\s*([\d.]+)\s*(mm|cm|in)?\s*", str(value).lower())
    if not m:
        raise ValueError(f"not a length: {value!r} (use e.g. 170mm, 17cm or 6.7in)")
    return float(m.group(1)) * {"mm": 1.0, "cm": 10.0, "in": 25.4}[m.group(2) or "mm"]


def parse_page(value) -> tuple[float, float]:
    """'letter', 'a4', 'nature-communications' or 'WxH' with an optional unit ('210x279mm', '8.5x11in') -> (w, h) in mm."""
    if isinstance(value, (tuple, list)):
        return float(value[0]), float(value[1])
    v = str(value).strip().lower()
    if v in PAGE_SIZES_MM:
        return PAGE_SIZES_MM[v]
    m = re.fullmatch(r"([\d.]+)\s*x\s*([\d.]+)\s*(mm|cm|in)?", v)
    if not m:
        raise ValueError(f"not a page size: {value!r} (use letter, a4, nature-communications or e.g. 210x279mm, 8.5x11in)")
    f = {"mm": 1.0, "cm": 10.0, "in": 25.4}[m.group(3) or "mm"]
    return float(m.group(1)) * f, float(m.group(2)) * f


def figure_label(style: JournalStyle, *, number=None, supplementary: bool = False, prefix: str | None = None) -> str:
    """'Supplementary Fig. 4' (nature-communications) or 'Supplementary Figure 4' (generic); ``prefix`` overrides the word."""
    if number in (None, ""):
        return ""
    word = prefix if prefix is not None else (style.supplementary_prefix if supplementary else style.figure_prefix)
    return f"{word} {number}".strip()


def panel_letters(n: int, case: str = "lower") -> list[str]:
    letters = list("abcdefghijklmnopqrstuvwxyz"[:n])
    return [x.upper() for x in letters] if case == "upper" else letters


def use_fonts(families=SANS_SERIF, font_dirs=()) -> str:
    """Register fonts from ``font_dirs`` and set matplotlib to the first available family (text embedded as TrueType)."""
    for d in font_dirs or ():
        for f in sorted(glob.glob(str(Path(d) / "*.[ot]tf")) + glob.glob(str(Path(d) / "*.ttc"))):
            try:
                fm.fontManager.addfont(f)
            except Exception:
                pass
    families = list(families)
    if "DejaVu Sans" not in families:
        families.append("DejaVu Sans")
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": families, "pdf.fonttype": 42, "ps.fonttype": 42,
                         "svg.fonttype": "none"})
    names = {f.name for f in fm.fontManager.ttflist}
    return next(f for f in families if f in names or f == "DejaVu Sans")


_MARKUP = re.compile(r"(\*\*\*.+?\*\*\*|\*\*.+?\*\*|\*[^*\s][^*]*?\*)")


def parse_markup(text: str) -> list[list[tuple[str, bool, bool]]]:
    """Split legend text into words; each word is a list of (text, bold, italic) runs written without spaces between them,
    so '*Genus species*.' keeps the period against the italic word."""
    units: list[list[tuple[str, bool, bool]]] = []
    glue = False
    for part in _MARKUP.split(text):
        if not part:
            continue
        bold = italic = False
        inner = part
        if part.startswith("***") and part.endswith("***") and len(part) > 6:
            bold = italic = True; inner = part[3:-3]
        elif part.startswith("**") and part.endswith("**") and len(part) > 4:
            bold = True; inner = part[2:-2]
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            italic = True; inner = part[1:-1]
        tokens = inner.split()
        if not tokens:
            glue = False
            continue
        for k, tok in enumerate(tokens):
            if k == 0 and glue and not inner[0].isspace() and units:
                units[-1].append((tok, bold, italic))
            else:
                units.append([(tok, bold, italic)])
        glue = not inner[-1].isspace()
    return units


def _layout_words(fig, units, width_in: float, size: float) -> list[list[tuple[float, str, bool, bool]]]:
    """Greedy line breaking with word widths measured by the renderer; returns lines of (x offset in inches, text, bold, italic)."""
    rend = fig.canvas.get_renderer()
    cache: dict = {}

    def width(text, bold, italic):
        key = (text, bold, italic)
        if key not in cache:
            t = fig.text(0, 0, text, fontsize=size, fontweight="bold" if bold else "normal", fontstyle="italic" if italic else "normal")
            cache[key] = t.get_window_extent(rend).width / fig.dpi
            t.remove()
        return cache[key]

    space = width("x x", False, False) - width("xx", False, False)
    lines, cur, x = [], [], 0.0
    for unit in units:
        runs, w = [], 0.0
        for text, bold, italic in unit:
            runs.append((w, text, bold, italic)); w += width(text, bold, italic)
        start = x + space if cur else 0.0
        if cur and start + w > width_in:
            lines.append(cur); cur, start = [], 0.0
        cur.extend((start + dx, t, b, i) for dx, t, b, i in runs)
        x = start + w
    if cur:
        lines.append(cur)
    return lines


def default_caption(letters: list[str], fields: dict, *, has_gff: bool, has_reads: bool, has_rnaseq: bool, has_variants: bool) -> str:
    """Legend body describing panels a-d with the numbers of this figure."""
    a, b, c, d = (f"**{x}**" for x in letters[:4])
    rings = []
    if has_rnaseq:
        rings.append("RNA-seq depth")
    if has_gff:
        rings.append("tRNA, protein-coding and rRNA genes")
    if has_variants:
        rings.append("variant columns (bar height, minor-allele or insertion fraction; red, majority differs from the assembly; "
                     "black, deletion; purple, insertion; orange, minor allele)")
    rings.append("AT content")
    if has_reads and "n_drawn" in fields:
        rings.append(f"and the {fields['n_drawn']} longest of {fields['n_reads']} long reads, as spirals (one ring per circle) or arcs "
                     "shaded 5' to 3', thicker at insertions and thinner at deletions of 10 bp or more, with dots at read bases that "
                     "differ from the assembly")
    text = f"{a} Circular map of the {fields['length']:,}-bp mitogenome (redwood), from the outside in: " + "; ".join(rings) + ". "
    text += f"{b} mtDNA content per NUMT locus, each nuclear base counted once; dashed lines, class thresholds. "
    text += f"{c} NUMT loci on the nuclear scaffolds"
    if fields.get("max_factor", 1.0) > 1.0:
        text += (f", drawn at least {fields['min_box_kb']:.0f} kb wide (true spans {fields['span_lo_kb']:.1f}-"
                 f"{fields['span_hi_kb']:.1f} kb)")
    text += "; height, class; color, identity (75% or less share the lowest color). "
    text += (f"{d} Mitogenome intervals of each locus at its identity, colored by class; dotted lines join alignments of one locus "
             "(separate pieces or placements on repeats)" + ("; annotation below" if has_gff else "") + ". "
             "NUMTs were detected in the nuclear assembly under the search thresholds used.")
    return text


# ------------------------------------------------------------------------------------------------------------- figure
def plot_journal_figure(*, mito_fasta: Path, loci: list[dict], chrom_lengths: dict[str, int], out_base: Path,
                        gff: Path | None = None, main_bam: Path | None = None, rnaseq_bam: Path | None = None,
                        variant_table: Path | None = None, style: JournalStyle | str = "nature-communications", layout: str = "page",
                        page=None, figure_width=None, legend_columns: int | None = None, panel_labels: str | None = None,
                        panel_label_pt: float | None = None, text_pt=None, label: str = "", title: str = "", species: str = "",
                        caption: str | None = None, caption_append: str = "", fields: dict | None = None, legend: bool = True,
                        dpi: int = 300, fileforms=("pdf", "png"), keep_figure: bool = False, **circular_kwargs) -> dict:
    """Draw the composite (a map with ring key, b mtDNA content per NUMT locus, c landscape, d catalog with the annotation).

    ``layout="figure"`` writes just the figure at the figure width; ``layout="page"`` places it on a page in the style's
    geometry with the legend underneath. Plot heights scale with the figure width; text stays at the style's sizes."""
    from .numts import _identity_colorbar, draw_numt_catalog, draw_numt_landscape, draw_numt_sizes
    from .renderer import draw_circular_plot, read_reference

    st = STYLES[style] if isinstance(style, str) else style
    if layout not in ("figure", "page"):
        raise ValueError("layout must be 'figure' or 'page'")
    reference = read_reference(Path(mito_fasta)); L = len(reference)
    Wmm = parse_length_mm(figure_width) if figure_width is not None else st.figure_width_mm
    W = Wmm * MM; s = Wmm / 170.0
    cols = int(legend_columns or st.legend_columns)
    letters = panel_letters(4, panel_labels or st.panel_label_case)
    label_pt = float(panel_label_pt or st.panel_label_pt)
    tpt = tuple(text_pt) if text_pt else st.text_pt
    longest = max(chrom_lengths.values()) if chrom_lengths else 1
    shown = {c for c, l in chrom_lengths.items() if l >= 0.02 * longest} | {r["chrom"] for r in loci}
    nrows = max(1, min(40, len(shown)))
    h_map, b_h = 3.30 * s, 2.70 * s
    h_land = max(0.75, 0.075 * nrows + 0.10)
    h_cat = 0.86 * s + 0.38
    top_block = max(0.02 + h_map + 0.06, 0.34 + b_h + 0.30) + 0.21
    fig_h = top_block + h_land + 0.56 + h_cat + 0.30
    if layout == "figure":
        pad = 2.0 * MM
        PW, PH = W + 2 * pad, fig_h + 2 * pad
        X0, top = pad, PH - pad
    else:
        pw, ph = parse_page(page) if page is not None else st.page_mm
        PW, PH = pw * MM, ph * MM
        X0, top = (PW - W) / 2, PH - st.top_mm * MM
    fig = plt.figure(figsize=(PW, PH), dpi=dpi)

    def ax_in(x, y_top, w, h):
        return fig.add_axes([x / PW, (y_top - h) / PH, w / PW, h / PH])

    b_top = top - 0.34
    ax_sz = ax_in(X0 + 0.715 * W, b_top, 0.285 * W, b_h)
    ax0 = ax_in(X0 + 0.010 * W, top - 0.02, 0.655 * W, h_map); ax0.set_anchor("NW")
    draw_circular_plot(ax0, length=L, reference=reference, gff=Path(gff) if gff else None, main_bam=Path(main_bam) if main_bam else None,
                       rnaseq_bam=Path(rnaseq_bam) if rnaseq_bam else None, variant_table=Path(variant_table) if variant_table else None,
                       text_pt=tpt, **circular_kwargs)
    for sp in ax0.spines.values():
        sp.set_visible(False)
    fig.canvas.draw()                                             # keep the map's outer labels inside the figure box
    rend = fig.canvas.get_renderer()
    boxes = [t.get_window_extent(rend) for t in ax0.texts if t.get_text()]
    if boxes:
        dx = max(0.0, X0 - min(b.x0 for b in boxes) / fig.dpi)
        dy = max(0.0, max(b.y1 for b in boxes) / fig.dpi - top)
        if dx or dy:
            pos = ax0.get_position(original=True)
            ax0.set_position([pos.x0 + dx / PW, pos.y0 - dy / PH, pos.width, pos.height])
    draw_numt_sizes(ax_sz, loci, L, title="")
    cursor = min(ax0.get_position(original=True).y0 * PH - 0.06, b_top - b_h - 0.30) - 0.21
    ax1 = ax_in(X0 + 0.105 * W, cursor, 0.82 * W, h_land)
    n = {k: sum(1 for r in loci if r["class"] == k) for k in ("full-length", "large", "fragment")}
    draw_numt_landscape(ax1, loci, chrom_lengths, label_fontsize=min(max(5.0, tpt[0]), tpt[1]),
                        title=f"NUMT landscape: {len(loci)} loci ({n['full-length']} full-length, {n['large']} large, "
                              f"{n['fragment']} fragments) on {nrows} sequences")
    _identity_colorbar(fig, ax1, fraction=0.02, label_size=min(6, tpt[1]), tick_size=max(5, tpt[0]))
    cursor -= h_land + 0.56
    ax2 = ax_in(X0 + 0.105 * W, cursor, 0.82 * W, h_cat)
    draw_numt_catalog(ax2, loci, L, gff, legend=False,
                      title=f"NUMT catalog: {len(loci)} loci, {sum(r.get('mtdna_bp', r['mito_bp']) for r in loci):,} bp of mtDNA in the nuclear genome")
    ax2.set_ylabel("identity (%)", fontsize=min(7, tpt[1]))
    figure_bottom = cursor - h_cat - 0.30
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    for letter, ax in zip(letters, (ax0, ax_sz, ax1, ax2)):
        box = ax.get_position(original=True)
        if ax is ax0:
            x, y, va = X0, top, "top"
        elif ax is ax_sz:
            x, y, va = box.x0 * PW - 0.50, top, "top"
        else:
            ttl = ax._left_title if ax._left_title.get_text() else ax.title
            x, y, va = X0, ttl.get_window_extent(rend).y0 / fig.dpi, "bottom"
        fig.text(x / PW, y / PH, letter, fontsize=label_pt, fontweight="bold", va=va, ha="left")

    counts = getattr(ax0, "_redwood_read_counts", None)
    data = {"length": L, **dict(zip("abcd", (f"**{x}**" for x in letters)))}
    data.update(getattr(ax1, "_redwood_landscape", {}))
    if counts:
        data.update({"n_reads": counts["reads"], "n_drawn": counts["multipass_drawn"] + counts["regular_drawn"]})
    data.update(fields or {})
    if not title and species:
        title = f"{DEFAULT_TITLE} of *{species}*."
    head = f"{label}{st.label_separator}{title}" if label and title else (label or title)
    body = caption.format(**data) if caption is not None else default_caption(
        letters, data, has_gff=bool(gff), has_reads=bool(main_bam), has_rnaseq=bool(rnaseq_bam), has_variants=bool(variant_table or main_bam))
    if caption_append:
        body = f"{body} {caption_append.format(**data)}".strip()
    markup = (f"**{head}** " if head else "") + body
    result = {"outputs": [], "legend_markdown": None, "warnings": [], "read_counts": counts,
              "layout_mm": {"page": [round(PW / MM, 1), round(PH / MM, 1)], "figure_left": round(X0 / MM, 1), "figure_width": round(Wmm, 1),
                            "figure_top": round((PH - top) / MM, 1), "figure_bottom": round((PH - figure_bottom) / MM, 1)}}
    if layout == "page" and legend and markup.strip():
        units = [[(t, True, i) for t, _b, i in u] for u in parse_markup(f"{head}")] if head else []
        units += parse_markup(body)
        measure = min(st.measure_mm, PW / MM - 20.0)
        col_w = (measure - (cols - 1) * st.gutter_mm) / cols * MM
        lines = _layout_words(fig, units, col_w, st.legend_pt)
        per = math.ceil(len(lines) / cols)
        lead = st.leading_pt / 72
        x_start = (PW - measure * MM) / 2
        first_baseline = figure_bottom - st.legend_gap_mm * MM - st.legend_pt / 72 * 0.78
        bottom = first_baseline
        for col in range(cols):
            for k, line in enumerate(lines[col * per:(col + 1) * per]):
                yb = first_baseline - k * lead
                for dxw, text, bold, italic in line:
                    fig.text((x_start + col * (col_w + st.gutter_mm * MM) + dxw) / PW, yb / PH, text, fontsize=st.legend_pt,
                             fontweight="bold" if bold else "normal", fontstyle="italic" if italic else "normal", va="baseline", ha="left")
                bottom = min(bottom, yb)
        result["layout_mm"].update({"legend_lines": len(lines), "legend_columns": cols, "legend_bottom": round((PH - bottom) / MM, 1)})
        if bottom < st.bottom_mm * MM:
            result["warnings"].append(f"legend ends {(st.bottom_mm * MM - bottom) / MM:.1f} mm below the {st.name} text block; "
                                      "shorten the legend or use a narrower figure")
    out_base = Path(out_base)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    for ff in fileforms:
        p = Path(f"{out_base}.{ff}"); fig.savefig(p, dpi=dpi); result["outputs"].append(str(p))
    if legend and markup.strip():
        md = Path(f"{out_base}.legend.md"); md.write_text(markup + "\n"); result["legend_markdown"] = str(md)
    if keep_figure:
        result["figure"] = fig
    else:
        plt.close(fig)
    return result


def run_composite(args) -> dict:
    from .numts import read_loci

    st = STYLES[args.style]
    font = use_fonts(tuple(args.font or ()) + st.fonts, args.font_dir or ())
    fai = Path(str(args.nuclear_fasta) + ".fai")
    if not fai.exists():
        import pysam
        pysam.faidx(str(args.nuclear_fasta))
    lengths = {f[0]: int(f[1]) for f in (l.split("\t") for l in fai.read_text().splitlines() if l.strip())}
    label = args.figure_label or figure_label(st, number=args.figure_number, supplementary=args.supplementary, prefix=args.label_prefix)
    fields = {}
    for kv in args.field or []:
        k, _, v = kv.partition("=")
        fields[k] = v
    res = plot_journal_figure(
        mito_fasta=args.mito_fasta, loci=read_loci(args.numt_loci), chrom_lengths=lengths, out_base=args.output_base, gff=args.gff,
        main_bam=args.long_read_bam, rnaseq_bam=args.rnaseq_bam, variant_table=args.variant_table, style=st, layout=args.layout,
        page=args.page, figure_width=args.figure_width, legend_columns=args.legend_columns, panel_labels=args.panel_labels,
        panel_label_pt=args.panel_label_size, text_pt=args.text_size, label=label, title=args.title or "", species=args.species or "",
        caption=Path(args.caption_file).read_text().strip() if args.caption_file else None, caption_append=args.caption_append or "",
        fields=fields, legend=not args.no_legend, dpi=args.dpi, fileforms=tuple(args.fileform))
    res["font"] = font
    print(json.dumps(res, indent=2, default=str))
    for w in res["warnings"]:
        print("warning:", w)
    return res
