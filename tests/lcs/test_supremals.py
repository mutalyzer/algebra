import pytest
from algebra import Variant
from algebra.lcs.supremals import supremal, supremal_sequence


@pytest.mark.parametrize("reference, variants, expected", [
    ("GTGTGTTTTTTTAACAGGGA", [Variant(8, 9, "")], Variant(5, 12, "TTTTTT")),
    ("ACTG", [Variant(0, 1, "A")], Variant(0, 0, "")),
    ("TGCATTAGGGCAAGGGTCTTCGACTTTCCACGAAAATCGCGTCGGTTTGAC", [Variant(24, 25, "")], Variant(24, 27, "TT")),
])
def test_supremal(reference, variants, expected):
    variant, *_ = supremal(reference, variants, offset=1)
    assert variant == expected


@pytest.mark.parametrize("reference, observed, expected", [
    ("TGCATTAGGGCAAGGGTCTTCGACTTTCCACGAAAATCGCGTCGGTTTGAC", "TGCATTAGGGCAAGGGTCTTCGACTTCCACGAAAATCGCGTCGGTTTGAC", Variant(24, 27, "TT")),
])
def test_supremal_sequence(reference, observed, expected):
    variant, *_ = supremal_sequence(reference, observed)
    assert variant == expected

