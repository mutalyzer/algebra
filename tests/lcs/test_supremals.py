import pytest
from algebra.lcs.supremals import supremal, supremal_sequence, trim
from algebra.variants import Variant, patch


TESTS = [
    # reference, variants, supremal
    ("GTGTGTTTTTTTAACAGGGA", [Variant(8, 9, "")], Variant(5, 12, "TTTTTT")),
    ("ACTG", [Variant(0, 1, "A")], Variant(0, 0, "")),
    ("TGCATTAGGGCAAGGGTCTTCGACTTTCCACGAAAATCGCGTCGGTTTGAC", [Variant(24, 25, "")], Variant(24, 27, "TT")),
    ("TGCATTAGGGCAAGGGTCTTCGACTTTCCACGAAAATCGC", [Variant(24, 25, "")], Variant(24, 27, "TT")),
    ("GGGTCTTCGACTTTCCACGAAAATCGC", [Variant(11, 12, "")], Variant(11, 14, "TT")),
    ("AAA", [Variant(0, 1, "")], Variant(0, 3, "AA")),
    ("A", [Variant(0, 1, "")], Variant(0, 1, "")),
]


@pytest.mark.parametrize("reference, variants, expected", TESTS)
def test_supremal(reference, variants, expected):
    variant, *_ = supremal(reference, variants)
    assert variant == expected


@pytest.mark.parametrize("reference, variants, expected", TESTS)
def test_supremal_sequence(reference, variants, expected):
    variant, *_ = supremal_sequence(reference, patch(reference, variants))
    assert variant == expected


@pytest.mark.parametrize("reference, observed, prefix_len, suffix_len", [
    ("", "", 0, 0),
    ("A", "A", 1, 0),
    ("AA", "A", 1, 0),
    ("AAA", "AA", 2, 0),
    ("A", "C", 0, 0),
    ("AAATAAA", "T", 0, 0),
    ("AAATAAA", "AAACAAA", 3, 3),
    ("AAATAAA", "AAATAAA", 7, 0),
])
def test_trim(reference, observed, prefix_len, suffix_len):
    assert trim(reference, observed) == (prefix_len, suffix_len)
