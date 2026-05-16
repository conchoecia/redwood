from redwood.cli import build_parser


def test_plot_parser_defaults():
    parser = build_parser()
    args = parser.parse_args(["plot", "--gff", "annotation.gff"])

    assert args.command == "plot"
    assert args.gff.endswith("annotation.gff")
    assert args.fileform == ["png"]
    assert args.doubled == []
    assert args.transparent is True
    assert args.verbose is False
    assert args.mito_fasta is None
    assert args.max_reads == 80


def test_end_to_end_parser_accepts_references_and_reads():
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--mito-fasta",
            "mito.fa",
            "--nuclear-fasta",
            "nuclear.fa",
            "--long-reads",
            "ont.fastq.gz",
            "--rnaseq-reads",
            "rna_1.fastq.gz",
            "rna_2.fastq.gz",
            "--outdir",
            "redwood-out",
        ]
    )

    assert args.command == "run"
    assert str(args.mito_fasta) == "mito.fa"
    assert str(args.nuclear_fasta) == "nuclear.fa"
    assert [str(path) for path in args.long_reads] == ["ont.fastq.gz"]
    assert [str(path) for path in args.rnaseq_reads] == ["rna_1.fastq.gz", "rna_2.fastq.gz"]
    assert args.long_read_depth == 100.0
    assert args.max_reads == 80


def test_plot_parser_accepts_extra_tracks():
    parser = build_parser()
    args = parser.parse_args(
        [
            "plot",
            "--gff",
            "annotation.gff",
            "--mito-fasta",
            "mito.fa",
            "--extra-track",
            "at",
            "metrics",
        ]
    )

    assert args.extra_tracks == ["at", "metrics"]
    assert args.mito_fasta.endswith("mito.fa")


def test_advanced_parser_accepts_mapping_steps():
    parser = build_parser()
    args = parser.parse_args(
        [
            "advanced",
            "map-long",
            "--mito-fasta",
            "mito.fa",
            "--long-reads",
            "ont.fastq.gz",
            "--outdir",
            "redwood-work",
        ]
    )

    assert args.command == "advanced"
    assert args.advanced_command == "map-long"
    assert str(args.mito_fasta) == "mito.fa"
