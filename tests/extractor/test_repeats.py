import pytest

from algebra.extractor import extract_repeats

@pytest.mark.parametrize(
    "sequence, expected",
    [
        ("AAB", [
            (0, 1, 2, 0),
        ]),
        ("AABAB", [
            (0, 1, 2, 0),
            (1, 2, 2, 0),
        ]),
        ("AABC", [
            (0, 1, 2, 0),
        ]),
        ("AABAAB", [
            (0, 1, 2, 0),
            (0, 3, 2, 0),
            (3, 1, 2, 0),
        ]),
        ("AABAABA", [
            (0, 1, 2, 0),
            (0, 3, 2, 1),
            (3, 1, 2, 0),
        ]),
        ("AABAABB", [
            (0, 1, 2, 0),
            (0, 3, 2, 0),
            (3, 1, 2, 0),
            (5, 1, 2, 0),
        ]),
        ("AABBAABB", [
            (0, 1, 2, 0),
            (0, 4, 2, 0),
            (2, 1, 2, 0),
            (4, 1, 2, 0),
            (6, 1, 2, 0),
        ]),
        ("CATCATACATACTACTAAAAA", [
            (0, 3, 2, 0),
            (3, 4, 2, 1),
            (9, 3, 2, 2),
            (16, 1, 5, 0),
        ]),
        ("AAABAAAB", [
            (0, 1, 3, 0),
            (0, 4, 2, 0),
            (4, 1, 3, 0),
        ]),

        ("CATCATACATACTACTAAAAA", [
            (0, 3, 2, 0),
            (3, 4, 2, 1),
            (9, 3, 2, 2),
            (16, 1, 5, 0),
        ]),
        ("ACCACCTCT", [
            (0, 3, 2, 0),
            (1, 1, 2, 0),
            (4, 1, 2, 0),
            (5, 2, 2, 0),
        ]),
        ("CCACACAC", [
            (0, 1, 2, 0),
            (1, 2, 3, 1),
        ]),
        ("CACATCACATCACATCACAT", [
            (0, 2, 2, 0),
            (0, 5, 4, 0),
            (5, 2, 2, 0),
            (10, 2, 2, 0),
            (15, 2, 2, 0),
        ]),
        ("AACAACAAAACAACA", [
            (0, 1, 2, 0),
            (0, 3, 2, 2),
            (3, 1, 2, 0),
            (3, 5, 2, 0),
            (6, 1, 4, 0),
            (8, 3, 2, 1),
            (11, 1, 2, 0),
        ]),
        ("AAAAAACAAAAAACAAAAAACACACAAAAAAA", [
            (0, 1, 6, 0),
            (0, 7, 3, 1),
            (7, 1, 6, 0),
            (14, 1, 6, 0),
            (19, 2, 3, 1),
            (25, 1, 7, 0),
        ]),
    ],


)
def test_repeats(sequence, expected):
    assert expected == extract_repeats(sequence)
