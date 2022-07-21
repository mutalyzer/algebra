import pytest
from algebra import Relation
from algebra.lcs import edit, lcs_graph
from algebra.supremal import (compare, find_supremal, spanning_variant,
                              subtract, union, intersect)
from algebra.variants import Variant, patch


@pytest.mark.parametrize("reference, lhs, rhs, expected", [
    ("AGATCCATTGTCAATGACAT", Variant(7, 11, "T"), Variant(10, 12, "CC"), Relation.OVERLAP),
    ("ATGCTATCCCTCCCCACTCC", Variant(7, 10, "CC"), Variant(10, 15, "TTCCC"), Relation.DISJOINT),
    ("GTGTGTTTTTTTAACAGGGA", Variant(6, 6, "G"), Variant(4, 12, "GGTTTTTTTT"), Relation.DISJOINT),
    ("CCACC", Variant(2, 3, "T"), Variant(2, 2, "T"), Relation.CONTAINS),
    ("CCACC", Variant(2, 3, "T"), Variant(3, 3, "T"), Relation.CONTAINS),
    ("CCACC", Variant(2, 2, "T"), Variant(2, 2, "T"), Relation.EQUIVALENT),
    ("CCACC", Variant(1, 1, "T"), Variant(2, 2, "T"), Relation.DISJOINT),
])
def test_compare(reference, lhs, rhs, expected):
    assert compare(reference, lhs, rhs) == expected


@pytest.mark.parametrize("reference, observed, expected_variant", [
    ("A", "C", Variant(0, 1, "C")),
    ("", "A", Variant(0, 0, "A")),
    ("A", "", Variant(0, 1)),
    ("ACCCA", "ACCA", Variant(1, 4, "CC")),
    ("AACCT", "AACGT", Variant(2, 4, "CG")),
])
def test_spanning_variant(reference, observed, expected_variant):
    _, lcs_nodes = edit(reference, observed)
    _, edges = lcs_graph(reference, observed, lcs_nodes)
    assert spanning_variant(reference, observed, edges) == expected_variant


@pytest.mark.parametrize("reference, observed, exception, message", [
    ("", "", ValueError, "No variants"),
])
def test_spanning_variant_fail(reference, observed, exception, message):
    _, lcs_nodes = edit(reference, observed)
    _, edges = lcs_graph(reference, observed, lcs_nodes)
    with pytest.raises(exception) as exc:
        spanning_variant(reference, observed, edges)
    assert str(exc.value) == message


@pytest.mark.parametrize("reference, variant, supremal_variant", [
    ("GTGTGTTTTTTTAACAGGGA", Variant(8, 9), Variant(5, 12, "TTTTTT")),
])
def test_find_supremal(reference, variant, supremal_variant):
    assert find_supremal(reference, variant, offset=1) == supremal_variant


@pytest.mark.parametrize("reference, lhs, rhs, expected", [
    ("TGGA", Variant(0, 3, "GGG"), Variant(1, 3, "GGG"), [Variant(0, 1, "")]),
    ("ACCTGC", Variant(1, 5, "CCCT"), Variant(3, 5, "TT"), [Variant(1, 2, "C"), Variant(2, 3, "C"), Variant(3, 4, "C")]),
    ("CGTCCA", Variant(1, 5, "TCCC"), Variant(1, 3, "TT"), [Variant(2, 3, "C"), Variant(3, 4, "C"), Variant(4, 5, "C")]),
    ("ACCTGC", Variant(1, 5, "CCCT"), Variant(1, 4, "CCC"), [Variant(4, 5, "T")]),
    ("ACTG", Variant(2, 4), Variant(2, 3), [Variant(3, 4)]),
    ("ACCTGC", Variant(1, 5, "TCTT"), Variant(3, 5, "TT"), [Variant(1, 2, "T"), Variant(2, 3, "C")]),
    ("ACCTGC", Variant(1, 5, "TCTT"), Variant(1, 3, "TC"), [Variant(3, 4, "T"), Variant(4, 5, "T")]),
    ("CACACAC", Variant(1, 6, "TCTCT"), Variant(3, 4, "T"), [Variant(1, 2, "T"), Variant(2, 3, "C"), Variant(4, 5, "C"), Variant(5, 6, "T")]),
    ("CACACAC", Variant(1, 6, "TCTCT"), Variant(1, 6, "TCACT"), [Variant(3, 4, "T")]),
    ("ACCACATTA", Variant(3, 4, "AGGGTA"), Variant(3, 4, "AGGGTA"), []),
    ("ACCACATTA", Variant(3, 4, "AGGGTA"), Variant(3, 4, "AGG"), [Variant(4, 4, "G"), Variant(4, 4, "T"), Variant(4, 4, "A")]),
    ("ACCACATTA", Variant(3, 4, "AGGGTA"), Variant(3, 4, "AGTA"), [Variant(4, 4, "G"), Variant(4, 4, "G")]),
    ("T", Variant(0, 1, "G"), Variant(0, 1), [Variant(1, 1, "G")]),
    ("TT", Variant(0, 2), Variant(0, 2, "T"), [Variant(0, 1)]),
])
def test_subtract(reference, lhs, rhs, expected):
    assert list(subtract(reference, lhs, rhs)) == expected


@pytest.mark.parametrize("reference, lhs, rhs", [
    ("CATATATC", Variant(1, 7, "ATATATAT"), Variant(4, 5, "TT")),
    ("AA", Variant(1, 1, "C"), Variant(1, 1, "T")),
    ("G", Variant(0, 1, "ATGG"), Variant(0, 0, "T")),
    ("G", Variant(0, 1, "TAGG"), Variant(0, 0, "T")),
    ("CCG", Variant(0, 3, "TAT"), Variant(1, 1, "A")),
    ("TT", Variant(0, 2, "T"), Variant(0, 2)),
])
def test_subtract_fail(reference, lhs, rhs):
    with pytest.raises(ValueError) as exc:
        print(list(subtract(reference, lhs, rhs)))
    assert str(exc.value) == "Undefined"


@pytest.mark.parametrize("reference, lhs, rhs, expected", [
])
def test_union(reference, lhs, rhs, expected):
    assert list(union(reference, lhs, rhs)) == expected


@pytest.mark.parametrize("reference, lhs, rhs", [
])
def test_union_fail(reference, lhs, rhs):
    with pytest.raises(ValueError) as exc:
        list(union(reference, lhs, rhs))
    assert str(exc.value) == "Undefined"


@pytest.mark.parametrize("reference, lhs, rhs, expected", [
])
def test_intersect(reference, lhs, rhs, expected):
    assert list(intersect(reference, lhs, rhs)) == expected
