from dataclasses import replace
import json

import numpy as np
import pysam
import pytest

from redwood.cli import build_parser, main
from redwood.linear import feature_polygon, rna_depth
from redwood.linear_evidence import (
    Reference, classify_junction, collect_evidence, load_annotations, load_classes,
    load_reference, read_segments, select_display_segments, write_evidence,
)
from redwood.workflow import map_long, prepare_reference, run_end_to_end


@pytest.fixture
def fixture(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">mito descriptive header\n" + "ACGT" * 250 + "\n")
    return tmp_path, fasta, Reference("mito", "ACGT" * 250)


def write_bam(directory, rows, name="reads", length=1000, contig="mito"):
    raw, path = directory / f"{name}.unsorted.bam", directory / f"{name}.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": contig, "LN": length}]}
    with pysam.AlignmentFile(raw, "wb", header=header) as bam:
        for name, start, cigar, flag in rows:
            read = pysam.AlignedSegment(bam.header)
            read.query_name = name
            read.reference_id = 0
            read.reference_start = start
            read.cigarstring = cigar
            read.flag = flag
            read.mapping_quality = 0 if name == "mapq0" else 60
            read.query_sequence = "A" * sum(n for op, n in read.cigartuples if op in {0, 1, 4, 7, 8})
            bam.write(read)
    pysam.sort("-o", str(path), str(raw))
    pysam.index(str(path))
    return path


def test_linear_coordinates_depth_clips_and_classes(fixture):
    directory, _, reference = fixture
    bam = write_bam(directory, [
        ("full", 0, "1000M", 0),
        ("gaps", 0, "5S3=2X2D2M3N2M5S", 16),
        ("full", 990, "990S10M", 2048),  # overlaps full's primary; one molecule in depth
        ("secondary", 0, "1000M", 256),
        ("mapq0", 999, "1M100S", 0),
    ])
    segments, counts = read_segments(bam, reference, {"full": {"class": "mito", "locus": ""}})
    evidence = collect_evidence(segments, reference, counts)
    np.testing.assert_array_equal(evidence["depth"]["mito"], np.ones(1000))
    other = evidence["depth"]["unclassified"]
    assert other[:15].tolist() == [1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0]
    assert other[999] == 1
    assert counts == {"primary": 3, "supplementary": 1, "secondary": 1, "hard_clipped": 0}
    summary = evidence["summary"]
    assert summary["reads"] == 3
    assert summary["terminal_starts"] == 2
    assert summary["terminal_ends"] == 3
    assert summary["terminal_soft_clips"]["left"]["short"] == 1
    assert summary["terminal_soft_clips"]["right"]["long"] == 1
    assert summary["mapq0_alignment_fraction"] == .25
    assert evidence["profiles"]["mapq0_starts"][-1] == 1
    assert summary["topology_call"] == "not_assessed"
    write_evidence(directory / "evidence", evidence)
    rows = (directory / "evidence.profiles.tsv").read_text().splitlines()
    assert rows[1].startswith("1\t") and rows[-1].startswith("1000\t")
    assert len(rows) == 1001


def test_split_junction_orientation_and_breakpoints(fixture):
    directory, _, reference = fixture
    bam = write_bam(directory, [
        ("forward", 900, "100M100S", 0), ("forward", 0, "100S100M", 2048),
        ("reverse", 0, "100S100M", 16), ("reverse", 900, "100M100S", 2064),
        ("fold", 900, "100M100S", 0), ("fold", 900, "100M100S", 2064),
        # Both alignments reach a terminus, but the query-facing breakpoints are internal.
        ("internal", 0, "500M500S", 0), ("internal", 0, "500M500S", 2064),
    ])
    segments, counts = read_segments(bam, reference)
    evidence = collect_evidence(segments, reference, counts, junction_window=100)
    kinds = {row["read"]: row["kind"] for row in evidence["junctions"]}
    assert kinds == {"forward": "end_to_start", "reverse": "end_to_start",
                     "fold": "terminal_inverted", "internal": "internal_inverted"}
    left = next(s for s in segments if s.name == "forward" and not s.supplementary)
    right = next(s for s in segments if s.name == "forward" and s.supplementary)
    assert classify_junction(left, replace(right, query_start=1000), 1000, 100) == "unassessed_query_gap_or_overlap"
    assert classify_junction(left, replace(right, query_start=-1000), 1000, 100) == "unassessed_query_gap_or_overlap"


def test_hard_clips_no_supplementaries_and_soft_masking(fixture):
    directory, _, reference = fixture
    bam = write_bam(directory, [("hard", 0, "300H1000M", 0)])
    reference = replace(reference, sequence=reference.sequence.lower())
    segments, counts = read_segments(bam, reference)
    evidence = collect_evidence(segments, reference, counts)
    summary = evidence["summary"]
    assert summary["lowercase_fraction"] == 1
    assert summary["reads_with_300bp_unaligned"] == 1
    assert summary["terminal_soft_clips"]["left"]["unassessed_hard_clip"] == 1
    assert "Hard-clipped" in " ".join(summary["warnings"])
    assert "No supplementary" in " ".join(summary["warnings"])
    assert classify_junction(segments[0], segments[0], 1000, 100) == "unassessed_hard_clip"


def test_input_identity_and_doubled_rejection(fixture):
    directory, fasta, reference = fixture
    bam = write_bam(directory, [], length=2000)
    with pytest.raises(ValueError, match="length mismatch"):
        read_segments(bam, reference)
    with pytest.raises(ValueError, match="does not contain"):
        read_segments(bam, replace(reference, name="wrong"))
    fasta.write_text(">one\nACGT\n>two\nACGT\n")
    with pytest.raises(ValueError, match="exactly one"):
        load_reference(fasta)
    with pytest.raises(ValueError, match="cannot be combined"):
        main(["plot", "--topology", "linear", "--doubled", "main"])


def test_annotations_deduplicate_and_clip_without_wrapping(fixture):
    directory, _, reference = fixture
    gff = directory / "ref.gff"
    gff.write_text("##sequence-region mito 1 1000\n"
                   "mito\ttest\tgene\t10\t200\t.\t+\t.\tName=parent\n"
                   "mito\ttest\tCDS\t10\t200\t.\t+\t0\tName=COX1\n"
                   "mito\ttest\tgene\t-10\t8\t.\t-\t.\tName=partial\n"
                   "mito\ttest\trepeat_region\t950\t1100\t.\t-\t.\tName=ITR\n"
                   "##FASTA\n>mito\nACGT\n")
    features = load_annotations(gff, reference)
    assert [f["name"] for f in features] == ["COX1", "partial", "ITR"]
    assert features[1]["start"] == -11 and features[2]["stop"] == 1100
    polygon = feature_polygon(949, 1100, 1, "-", 1000)
    assert min(x for x, _ in polygon) == 949.5
    assert max(x for x, _ in polygon) == 1000.5
    with pytest.raises(ValueError, match="sequence-region"):
        load_annotations(gff, replace(reference, sequence="ACGT"))
    gff.write_text("other\ttest\tgene\t1\t200\t.\t+\t.\tName=bad\n")
    with pytest.raises(ValueError, match="does not match"):
        load_annotations(gff, reference)


def test_query_and_display_limit_do_not_change_evidence(fixture):
    directory, _, reference = fixture
    bam = write_bam(directory, [("full", 0, "1000M", 0), ("clipped", 0, "1000M50S", 0)])
    segments, counts = read_segments(bam, reference)
    selected = select_display_segments(segments, 1000, ["ALNLEN >= 1000", "MAPLEN <= reflength"])
    assert [s.name for s in selected] == ["full"]
    assert len(select_display_segments(segments, 1000)) == 2  # no legacy default excludes full-length molecules
    assert select_display_segments(segments, 1000, max_reads=0) == []
    assert collect_evidence(segments, reference, counts)["summary"]["reads"] == 2
    with pytest.raises(ValueError, match="nonnegative"):
        select_display_segments(segments, 1000, max_reads=-1)


def test_classes_accept_original_phase_table_and_validate(fixture):
    directory, _, reference = fixture
    path = directory / "classes.tsv"
    path.write_text("read\tclass\treadlen\tnuclear_locus\nread1\tNUMT\t1000\tchr2:1-1000\n")
    assert load_classes(path) == {"read1": {"class": "NUMT", "locus": "chr2:1-1000"}}
    path.write_text("read\tclass\nread1\tmito\nread1\tNUMT\n")
    with pytest.raises(ValueError, match="conflicting"):
        load_classes(path)
    path.write_text("read\twrong\n")
    with pytest.raises(ValueError, match="header"):
        load_classes(path)
    path.write_text("read\tclass\n")
    bam = write_bam(directory, [("missing", 0, "1000M", 0)])
    segments, _ = read_segments(bam, reference, load_classes(path))
    assert segments[0].read_class == "unclassified"


@pytest.mark.parametrize("style", ["redwood", "diagnostic"])
def test_linear_renderer_and_rna_never_wrap(fixture, style):
    directory, fasta, reference = fixture
    bam = write_bam(directory, [("full", 0, "1000M", 0), ("reverse", 998, "1M1N", 16)])
    forward, reverse = rna_depth(bam, reference)
    assert forward[0] == 1 and reverse[998] == 1 and reverse[999] == 0 and reverse[0] == 0
    assert main(["plot", "--topology", "linear", "--mito-fasta", str(fasta),
                 "--linear-style", style,
                 "--main-bam", str(bam), "--rnaseq-bam", str(bam), "--depth-scale", "log",
                 "--fileform", "png", "svg", "pdf", "--no-timestamp", "--dpi", "60",
                 "-o", str(directory / "linear")]) == 0
    for ext in ("png", "pdf", "svg", "evidence.json", "junctions.tsv", "profiles.tsv"):
        assert (directory / f"linear.{ext}").stat().st_size > 0
    summary = json.loads((directory / "linear.evidence.json").read_text())
    assert summary["display"]["reads"] == 2
    assert summary["display"]["style"] == style


def test_production_composition_windows_do_not_join_termini():
    from redwood.linear_redwood import linear_base_fraction

    profile = linear_base_fraction("A" * 101 + "G" * 101)
    assert len(profile) == 202
    assert profile[0] == 1 and profile[-1] == 0
    assert profile[100] == pytest.approx(101 / 201)
    np.testing.assert_array_equal(linear_base_fraction("ATGC"), [.5] * 4)
    np.testing.assert_array_equal(linear_base_fraction("GGCC", "GC"), [1] * 4)


def test_production_packing_preserves_split_segments_and_read_indels(fixture):
    from redwood.linear_redwood import pack_linear_reads, read_outline

    directory, _, reference = fixture
    bam = write_bam(directory, [("split", 0, "100M100S", 0),
                                ("split", 900, "100S100M", 2048),
                                ("middle", 400, "100M", 0),
                                ("overlap", 980, "5M20I5M", 0)])
    segments, _ = read_segments(bam, reference)
    placed, rows = pack_linear_reads(segments, reference.length)
    assert rows == 2
    assert {name: row for name, _, row in placed} == {"split": 0, "middle": 0, "overlap": 1}
    assert len(next(group for name, group, _ in placed if name == "split")) == 2
    segment = next(s for s in segments if s.name == "overlap")
    outline = read_outline(segment, 1, .02, 10)
    assert outline[:, 0].min() == 980.5
    assert outline[:, 0].max() == 990.5


def test_production_gradient_preserves_indels_without_slicing_reads(fixture):
    import io
    import xml.etree.ElementTree as ET
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    from redwood.linear_redwood import add_gradient_read, read_outline
    from redwood.renderer import REDWOOD_GRADIENT

    directory, _, reference = fixture
    bam = write_bam(directory, [("indels", 100, "30M20I30M20D20M", 0)])
    segments, _ = read_segments(bam, reference)
    segment = segments[0]
    path = Path(read_outline(segment, 1, .02, 10))
    # The existing insertion expansion and deletion narrowing survive clipping.
    assert path.contains_point((130, 1.012))
    assert not path.contains_point((110, 1.012))
    assert path.contains_point((170, 1))
    assert not path.contains_point((170, 1.006))
    assert path.contains_point((190, 1.006))
    fig, ax = plt.subplots()
    ax.set(xlim=(100, 201), ylim=(.97, 1.03))
    mesh = add_gradient_read(ax, segment, 1, .02, 10, REDWOOD_GRADIENT)
    mesh.set_gid("gradient_read")
    buffer = io.StringIO()
    fig.savefig(buffer, format="svg")
    plt.close(fig)
    root = ET.fromstring(buffer.getvalue())
    ns = {"svg": "http://www.w3.org/2000/svg"}
    group = root.find(".//svg:g[@id='gradient_read']", ns)
    assert len(group.findall(".//svg:linearGradient", ns)) == 1
    assert not group.findall(".//svg:image", ns)
    assert len(group.findall(".//svg:path", ns)) == 1


def test_production_optional_bands_and_class_colors(fixture):
    directory, fasta, _ = fixture
    bam = write_bam(directory, [("full", 0, "1000M", 0)])
    classes = directory / "classes.tsv"
    classes.write_text("read\tclass\nfull\tmito\n")
    assert main(["plot", "--topology", "linear", "--mito-fasta", str(fasta),
                 "--main-bam", str(bam), "--read-classes", str(classes), "--read-color", "class",
                 "--linear-track", "depth", "--linear-track", "ends", "--linear-track", "clips",
                 "--show-terminal-sequences", "--extra-track", "gc", "--dark",
                 "--no-timestamp", "--dpi", "60", "-o", str(directory / "styled")]) == 0
    summary = json.loads((directory / "styled.evidence.json").read_text())
    assert summary["production_style"]["tracks"] == ["depth", "ends", "clips"]
    assert summary["production_style"]["read_color"] == "class"
    assert summary["production_style"]["read_rows"] == 1


def test_annotation_only_and_empty_bam_are_supported(fixture):
    directory, fasta, _ = fixture
    assert main(["plot", "--topology", "linear", "--mito-fasta", str(fasta),
                 "--no-timestamp", "--dpi", "60", "-o", str(directory / "reference")]) == 0
    bam = write_bam(directory, [])
    assert main(["plot", "--topology", "linear", "--mito-fasta", str(fasta),
                 "--main-bam", str(bam), "--no-timestamp", "--dpi", "60", "--hide-evidence",
                 "-o", str(directory / "empty")]) == 0
    summary = json.loads((directory / "empty.evidence.json").read_text())
    assert summary["reads"] == 0 and summary["read_length_p99"] is None


def test_linear_workflow_preserves_segments_and_single_reference(fixture, monkeypatch):
    directory, fasta, reference = fixture
    commands = []
    monkeypatch.setattr("redwood.workflow.run_pipeline", lambda cmds, dry_run: commands.extend(cmds))
    args = build_parser().parse_args(["advanced", "map-long", "--topology", "linear",
                                     "--mito-fasta", str(fasta), "--long-reads", "unused.fastq",
                                     "--outdir", str(directory / "work"), "--dry-run"])
    result = map_long(args)
    assert result["copies"] == 1
    assert load_reference(directory / "work/references/mitochondrion.linear.fa") == reference
    assert "-Y" in commands[0] and "--secondary=no" in commands[0]
    args.copies = "2"
    with pytest.raises(ValueError, match="one reference copy"):
        map_long(args)
    prepare = build_parser().parse_args(["advanced", "prepare-reference", "--topology", "linear",
                                        "--mito-fasta", str(fasta), "--outdir", str(directory / "prepared")])
    prepare_reference(prepare)
    assert load_reference(directory / "prepared/mitochondrion.linear.fa") == reference
    assert not (directory / "prepared/mitochondrion.doubled.fa").exists()


def test_linear_run_passes_topology_to_plot_and_metrics(fixture, monkeypatch):
    directory, fasta, _ = fixture
    bam = write_bam(directory, [("full", 0, "1000M", 0)])
    def mapped(args):
        assert args.topology == "linear"
        return {"output_bam": str(bam)}
    monkeypatch.setattr("redwood.workflow.map_long", mapped)
    args = build_parser().parse_args(["run", "--topology", "linear", "--mito-fasta", str(fasta),
                                     "--long-reads", "unused.fastq", "--outdir", str(directory / "run"),
                                     "--dpi", "60", "--max-reads", "0"])
    result = run_end_to_end(args)
    assert result["metrics"]["topology"] == "linear"
    summary = json.loads((directory / "run/redwood.evidence.json").read_text())
    assert summary["reads"] == 1 and summary["display"]["reads"] == 0
    assert (directory / "run/redwood.png").exists()
