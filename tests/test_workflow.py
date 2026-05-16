from redwood.workflow import add_circular_interval_depth


def test_add_circular_interval_depth_wraps_once():
    depth = [0] * 10

    add_circular_interval_depth(depth, 8, 13, 10)

    assert depth == [1, 1, 1, 0, 0, 0, 0, 0, 1, 1]


def test_add_circular_interval_depth_caps_full_spans():
    depth = [0] * 10

    add_circular_interval_depth(depth, 8, 23, 10)

    assert depth == [1] * 10
