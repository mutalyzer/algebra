import pytest
from algebra import LCSgraph, Relation, Variant, compare


def test_compare():
    reference = "AAAAA"
    lhs = LCSgraph.from_variants(reference, [Variant(0, 0, "T")])
    rhs = LCSgraph.from_variants(reference, [Variant(0, 0, "TA")])
    assert compare(reference, lhs, rhs) == Relation.is_contained
