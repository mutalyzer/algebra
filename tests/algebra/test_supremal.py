import pytest
from algebra import Relation
from algebra.lcs import edit, lcs_graph
from algebra.supremal import (compare, find_supremal, spanning_variant,
                              subtract, union)
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


@pytest.mark.parametrize("reference, lhs, rhs_a, rhs_b", [
    ("CACACAC", Variant(1, 6, "TCTCT"), Variant(1, 6, "TCACT"), Variant(3, 4, "T")),
    ("ACCTGC", Variant(1, 5, "TCTT"), Variant(1, 3, "TC"), Variant(3, 5, "TT")),
])
def test_subtract(reference, lhs, rhs_a, rhs_b):
    assert patch(reference, subtract(reference, lhs, rhs_a)) == patch(reference, [rhs_b])
    assert patch(reference, subtract(reference, lhs, rhs_b)) == patch(reference, [rhs_a])


@pytest.mark.parametrize("reference, lhs, rhs", [
    ("CATATATC", Variant(1, 7, "ATATATAT"), Variant(5, 6, "AA")),
    ("CATATATC", Variant(1, 7, "ATATATAT"), Variant(4, 5, "TT")),
    ("ACCTGC", Variant(1, 3, "TC"), Variant(3, 5, "TT")),
    ("CAC", Variant(1, 2, "T"), Variant(1, 2, "G")),
])
def test_subtract_fail(reference, lhs, rhs):
    with pytest.raises(ValueError) as exc:
        subtract(reference, lhs, rhs)
    assert str(exc.value) == "Undefined"


@pytest.mark.parametrize("reference, a, b, c", [
    ("CACACAC", Variant(1, 6, "TCTCT"), Variant(1, 6, "TCACT"), Variant(3, 4, "T")),
])
def test_union(reference, a, b, c):
    assert patch(reference, union(reference, b, c)) == patch(reference, [a])
    assert patch(reference, union(reference, c, b)) == patch(reference, [a])
    assert patch(reference, union(reference, a, b)) == patch(reference, [a])
    assert patch(reference, union(reference, a, c)) == patch(reference, [a])
    assert patch(reference, union(reference, a, a)) == patch(reference, [a])


@pytest.mark.parametrize("reference, a, b", [
    ("CAC", Variant(1, 2, "T"), Variant(1, 2, "G")),
])
def test_union_fail(reference, a, b):
    with pytest.raises(ValueError) as exc:
        union(reference, a, b)
    assert str(exc.value) == "Undefined"
