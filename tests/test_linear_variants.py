from types import SimpleNamespace

import matplotlib.pyplot as plt
import pysam
import pytest

from redwood.linear_evidence import Reference, read_segments
from redwood.linear_variants import collect_variant_calls, draw_variant_details, load_variant_sites


def variant_bam(tmp_path, rows, length=100):
    raw, path = tmp_path / "raw.bam", tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "mito", "LN": length}]}
    with pysam.AlignmentFile(raw, "wb", header=header) as bam:
        for name, start, cigar, sequence, quality, flag in rows:
            read = pysam.AlignedSegment(bam.header)
            read.query_name, read.reference_id, read.reference_start = name, 0, start
            read.cigarstring, read.flag, read.mapping_quality = cigar, flag, 60
            if sequence is not None:
                read.query_sequence = sequence
            if quality is not None:
                read.query_qualities = pysam.qualitystring_to_array(quality)
            bam.write(read)
    pysam.sort("-o", str(path), str(raw))
    pysam.index(str(path))
    return path


def test_variant_calls_deduplicate_conflicts_quality_and_reference_orientation(tmp_path):
    reference = Reference("mito", "ACGT" * 25)
    rows = [
        ("duplicate", 0, "10M", "A" * 10, "I" * 10, 0),
        ("duplicate", 0, "10M", "A" * 10, "I" * 10, 2048),
        ("conflict", 0, "10M", "A" * 10, "I" * 10, 0),
        ("conflict", 0, "10M", "C" * 10, "I" * 10, 2048),
        ("low", 0, "10M", "G" * 10, "!" * 10, 0),
        ("reverse", 0, "10M", "T" * 10, "I" * 10, 16),
        ("deletion", 0, "3M2D5M", "G" * 8, "I" * 8, 0),
        ("skip", 0, "3M2N5M", "G" * 8, "I" * 8, 0),
        ("missing", 0, "10M", None, None, 0),
        ("secondary", 0, "10M", "C" * 10, "I" * 10, 256),
        ("outside", 20, "10M", "C" * 10, "I" * 10, 0),
        ("fill", 0, "10M", "A" * 10, "!" * 10, 0),
        ("fill", 0, "10M", "C" * 10, "I" * 10, 2048),
        ("ambiguous", 0, "10M", "N" * 10, "I" * 10, 0),
        ("no_quality", 0, "10M", "A" * 10, None, 0),
    ]
    bam = variant_bam(tmp_path, rows)
    sites = [{"position": p, "reference": reference.sequence[p - 1], "label": ""} for p in (2, 5, 8)]
    calls, metadata = collect_variant_calls(bam, reference, sites)
    assert calls["duplicate"] == ["A"] * 3
    assert calls["conflict"] == ["conflict"] * 3
    assert calls["reverse"] == ["T"] * 3  # BAM sequence already faces the reference.
    assert calls["low"] == calls["no_quality"] == ["low_quality"] * 3
    assert calls["deletion"] == ["G", "deletion", "G"]
    assert calls["skip"] == ["G", "skip", "G"]
    assert calls["missing"] == ["missing_sequence"] * 3
    assert calls["outside"] == ["uncovered"] * 3
    assert calls["fill"] == ["C"] * 3
    assert calls["ambiguous"] == ["ambiguous"] * 3
    assert "secondary" not in calls
    assert metadata["population_reads"] == 11
    assert metadata["alignment_counts"]["secondary_excluded"] == 1
    site = metadata["sites"][1]
    assert site["called_reads"] == 3
    assert site["counts"]["A"] == site["counts"]["C"] == site["counts"]["T"] == 1
    assert site["allele_fractions"]["A"] == pytest.approx(1 / 3)
    assert site["counts"]["conflict"] == site["counts"]["deletion"] == site["counts"]["skip"] == 1


def test_variant_site_input_validation_and_source_table(tmp_path):
    reference = Reference("mito", "ACGT" * 25)
    path = tmp_path / "sites.tsv"
    path.write_text("pos\tv5_allele(hap1)\tminor_allele\tgene\tHiFi_minor_frac\n5\tA\tG\tCOX1\t0.99\n2\tC\tT\tND1\t0.12\n")
    assert load_variant_sites(path, reference) == [
        {"position": 2, "reference": "C", "label": "ND1"},
        {"position": 5, "reference": "A", "label": "COX1"},
    ]
    for contents, message in [("0\n", "outside"), ("101\n", "outside"), ("2\n2\n", "duplicate"),
                              ("position\tref\n2\tA\n", "does not match"), ("bad\n", "integer"),
                              ("position\n", "no sites")]:
        path.write_text(contents)
        with pytest.raises(ValueError, match=message):
            load_variant_sites(path, reference)
    path.write_text("# selected sites\n5\tlabel\n2\n")
    assert [site["position"] for site in load_variant_sites(path, reference)] == [2, 5]


def test_variant_companion_balances_capped_rows_but_counts_all_reads(tmp_path):
    reference = Reference("mito", "A" * 100)
    rows = [(f"left{index:02}", 0, "90M", "A" * 90, "I" * 90, 0) for index in range(20)]
    rows += [(f"right{index:02}", 10, "90M", "C" * 90, "I" * 90, 0) for index in range(20)]
    rows += [("unselected", 0, "100M", "T" * 100, "I" * 100, 0)]
    bam = variant_bam(tmp_path, rows)
    selected = [segment for segment in read_segments(bam, reference)[0] if segment.name != "unselected"]
    sites = tmp_path / "sites.tsv"
    sites.write_text("position\tlabel\n20\tone\n80\ttwo\n")
    args = SimpleNamespace(variant_sites=sites, main_bam=bam, width=3.5,
                           variant_min_base_quality=20, terminal_window=5,
                           linear_read_selection="terminal-balanced")
    fig = draw_variant_details(args, reference, selected)
    metadata = fig._redwood_variant_metadata
    assert metadata["population_reads"] == 41
    assert metadata["selected_count"] == 40
    assert metadata["displayed_count"] == 30
    assert sum(row["read"].startswith("left") for row in metadata["displayed_reads"]) == 15
    assert sum(row["read"].startswith("right") for row in metadata["displayed_reads"]) == 15
    assert metadata["sites"][0]["counts"]["T"] == 1
    assert metadata["sites"][0]["allele_fractions"]["A"] == pytest.approx(20 / 41)
    assert all(text.get_fontsize() >= 6 for text in fig.axes[0].texts)
    plt.close(fig)
    args.width = .5
    with pytest.raises(ValueError, match="increase --width"):
        draw_variant_details(args, reference, selected)
