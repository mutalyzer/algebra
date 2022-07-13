import pytest
from algebra import Relation, compare
from algebra.variants import Variant, patch


TESTS = [
    ("AGATCCATTGTCAATGACAT", Variant(7, 11, "T"), Variant(10, 12, "CC"), Relation.OVERLAP),
    ("ATGCTATCCCTCCCCACTCC", Variant(7, 10, "CC"), Variant(10, 15, "TTCCC"), Relation.DISJOINT),
    ("GTGTGTTTTTTTAACAGGGA", Variant(6, 6, "G"), Variant(4, 12, "GGTTTTTTTT"), Relation.DISJOINT),
    ("CCACC", Variant(2, 3, "T"), Variant(2, 2, "T"), Relation.CONTAINS),
    ("CCACC", Variant(2, 3, "T"), Variant(3, 3, "T"), Relation.CONTAINS),
    ("CCACC", Variant(2, 2, "T"), Variant(2, 2, "T"), Relation.EQUIVALENT),
]


@pytest.mark.parametrize("reference, lhs, rhs, expected", TESTS)
def test_compare(reference, lhs, rhs, expected):
    assert compare(reference, patch(reference, [lhs]), patch(reference, [rhs])) == expected
