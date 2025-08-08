import pytest
from algebra import LCSgraph, Relation, Variant, compare


@pytest.mark.parametrize("reference, lhs, rhs, expected", [
    ("CC", [Variant(0, 0, "C"), Variant(1, 2, "")], [Variant(0, 0, "C"), Variant(1, 1, "T")], Relation.disjoint),
    ("CTCCG", [Variant(4, 4, "TCC")], [Variant(0, 0, "GT")], Relation.overlap),
    ("CCCCC", [Variant(4, 5, "")], [Variant(0, 1, "")], Relation.equivalent),
    ("AACCCCTTGTTC", [Variant(2, 2, "C")], [Variant(6, 7, "C")], Relation.is_contained),
    ("AACATTTC", [Variant(1, 1, "AC")], [Variant(4, 5, "C")], Relation.overlap),
    ("AAA", [Variant(1, 1, "T")], [Variant(1, 1, "C")], Relation.disjoint),
])
def test_compare(reference, lhs, rhs, expected):
    lhs_graph = LCSgraph.from_variants(reference, lhs)
    rhs_graph = LCSgraph.from_variants(reference, rhs)
    assert compare(reference, lhs_graph, rhs_graph) == expected
