import pytest
from algebra import Relation, Variant
from algebra.lcs import edit, lcs_graph
from algebra.relations.supremal_based import (are_disjoint, are_equivalent, compare, contains, find_supremal,
                                              have_overlap, is_contained, spanning_variant)


TESTS = [
    ("AGATCCATTGTCAATGACAT", Variant(7, 11, "T"), Variant(10, 12, "CC"), Relation.OVERLAP),
    ("ATGCTATCCCTCCCCACTCC", Variant(7, 10, "CC"), Variant(10, 15, "TTCCC"), Relation.DISJOINT),
    ("GTGTGTTTTTTTAACAGGGA", Variant(6, 6, "G"), Variant(4, 12, "GGTTTTTTTT"), Relation.DISJOINT),
    ("CCACC", Variant(2, 3, "T"), Variant(2, 2, "T"), Relation.CONTAINS),
    ("CCACC", Variant(2, 3, "T"), Variant(3, 3, "T"), Relation.CONTAINS),
    ("CCACC", Variant(2, 2, "T"), Variant(2, 2, "T"), Relation.EQUIVALENT),
    ("CCACC", Variant(1, 1, "T"), Variant(2, 2, "T"), Relation.DISJOINT),
]


@pytest.mark.parametrize("reference, lhs, rhs, expected", TESTS)
def test_are_equivalent(reference, lhs, rhs, expected):
    assert are_equivalent(reference, lhs, rhs) == (expected == Relation.EQUIVALENT)


@pytest.mark.parametrize("reference, lhs, rhs, expected", TESTS)
def test_contains(reference, lhs, rhs, expected):
    assert contains(reference, lhs, rhs) == (expected == Relation.CONTAINS)


@pytest.mark.parametrize("reference, lhs, rhs, expected", TESTS)
def test_is_contained(reference, lhs, rhs, expected):
    assert is_contained(reference, lhs, rhs) == (expected == Relation.IS_CONTAINED)


@pytest.mark.parametrize("reference, lhs, rhs, expected", TESTS)
def test_have_overlap(reference, lhs, rhs, expected):
    assert have_overlap(reference, lhs, rhs) == (expected == Relation.OVERLAP)


@pytest.mark.parametrize("reference, lhs, rhs, expected", TESTS)
def test_are_disjoint(reference, lhs, rhs, expected):
    assert are_disjoint(reference, lhs, rhs) == (expected == Relation.DISJOINT)


@pytest.mark.parametrize("reference, lhs, rhs, expected", TESTS)
def test_compare(reference, lhs, rhs, expected):
    assert compare(reference, lhs, rhs) == expected


@pytest.mark.parametrize("reference, observed, expected", [
    ("A", "C", Variant(0, 1, "C")),
    ("", "A", Variant(0, 0, "A")),
    ("A", "", Variant(0, 1, "")),
    ("ACCCA", "ACCA", Variant(1, 4, "CC")),
    ("AACCT", "AACGT", Variant(2, 4, "CG")),
    ("", "", Variant(0, 0, "")),
])
def test_spanning_variant(reference, observed, expected):
    _, lcs_nodes = edit(reference, observed)
    _, edges = lcs_graph(reference, observed, lcs_nodes)
    assert spanning_variant(reference, observed, edges) == expected


@pytest.mark.parametrize("reference, variant, expected", [
    ("GTGTGTTTTTTTAACAGGGA", Variant(8, 9, ""), Variant(5, 12, "TTTTTT")),
    ("ACTG", Variant(0, 1, "A"), Variant(0, 0, "")),
])
def test_find_supremal(reference, variant, expected):
    supremal, *_ = find_supremal(reference, variant, offset=1)
    assert supremal == expected
