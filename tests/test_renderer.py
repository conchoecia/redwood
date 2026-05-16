from pathlib import Path

from redwood.renderer import plot_file


def test_renderer_writes_plot_from_example_dataset(tmp_path):
    repo = Path(__file__).resolve().parents[1]
    dataset = repo / "examples" / "datasets" / "human"
    output_base = tmp_path / "human"

    plot_file(
        output_base=str(output_base),
        fileforms=["png"],
        dpi=80,
        reference_fasta=dataset / "reference.fa",
        gff=dataset / "annotation.gff",
        main_bam=dataset / "reads.mapped.bam",
        rnaseq_bam=dataset / "rnaseq.mapped.bam",
        no_timestamp=True,
    )

    output = tmp_path / "human.png"
    assert output.exists()
    assert output.stat().st_size > 0
