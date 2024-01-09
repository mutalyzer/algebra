import pytest
from algebra.lcs.supremals import lcs_graph, lcs_graph_sequence, lcs_graph_supremal, trim
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
    ("A", [], Variant(0, 0, "")),
]


@pytest.mark.parametrize("reference, variants, expected", TESTS)
def test_lcs_graph(reference, variants, expected):
    graph = lcs_graph(reference, variants)
    assert graph.supremal == expected


@pytest.mark.parametrize("reference, variants, offset, expected", [
    ("XXXXXXXXXXCATATATCGXXXXXXXXXX", [Variant(11, 12, "T"), Variant(16, 17, "G"), Variant(18, 19, "AT")], 2, Variant(11, 19, "TTATAGCAT")),
    ("XXXXXXXXXXCATATATCGXXXXXXXXXX", [Variant(11, 12, "T"), Variant(16, 17, "G"), Variant(18, 19, "AT")], 3, Variant(11, 19, "TTATAGCAT")),
    ("XXXXXXXXXXCATATATCGXXXXXXXXXX", [Variant(11, 12, "T"), Variant(16, 17, "G"), Variant(18, 19, "AT")], 4, Variant(11, 19, "TTATAGCAT")),
    ("XXXXXXXXXXCATATATCGXXXXXXXXXX", [Variant(11, 12, "T"), Variant(16, 17, "G"), Variant(18, 19, "AT")], 40, Variant(11, 19, "TTATAGCAT")),
    ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", [Variant(20, 21, "T")], 1, Variant(0, 74, "AAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")),
])
def test_lcs_graph_offset(reference, variants, offset, expected):
    graph = lcs_graph(reference, variants, offset)
    assert graph.supremal == expected


@pytest.mark.parametrize("reference, variants, expected", TESTS)
def test_lcs_graph_sequence(reference, variants, expected):
    graph = lcs_graph_sequence(reference, patch(reference, variants))
    assert graph.supremal == expected


@pytest.mark.parametrize("reference, supremal", [(reference, supremal) for reference, _, supremal in TESTS])
def test_lcs_graph_supremal(reference, supremal):
    graph = lcs_graph_supremal(reference, supremal)
    assert graph.supremal == supremal


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
