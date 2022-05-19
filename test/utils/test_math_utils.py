from sp.utils import math_utils


def test_sum_numbers():
    test_cases = [
        ((3, 4), 7),
        ((1, 0), 1),
        ((0, 1), 1),
        ((-1, 1), 0),
    ]

    for (x, y), expected in test_cases:
        assert math_utils.add(x, y) == expected
