import pytest
from algebra import Relation, Variant
from algebra.lcs import edit, lcs_graph
from algebra.relations.supremal_based import compare, find_supremal, spanning_variant


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
