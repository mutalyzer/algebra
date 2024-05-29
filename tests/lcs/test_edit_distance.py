import pytest
from algebra.lcs import edit_distance


@pytest.mark.parametrize("reference, observed, expected", [
    ("CTCGGCATTA", "GGCTGGCTGT", 6),
    ("TTT", "TTTTAT", 3),
    ("TTTTAT", "TTT", 3),
])
def test_edit_distance(reference, observed, expected):
    assert edit_distance(reference, observed) == expected
