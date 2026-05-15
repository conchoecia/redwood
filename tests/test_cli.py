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
