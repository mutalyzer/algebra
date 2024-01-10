import pytest
from algebra.lcs import edit_distance


@pytest.mark.parametrize("reference, observed, expected", [
    ("CTCGGCATTA", "GGCTGGCTGT", 6),
])
def test_edit_distance(reference, observed, expected):
    assert edit_distance(reference, observed) == expected
