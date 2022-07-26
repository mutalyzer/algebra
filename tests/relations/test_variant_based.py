import pytest
from algebra import (Relation, Variant, are_disjoint, are_equivalent,
                     compare, contains, have_overlap, is_contained)


TESTS = [
    ("ACGT", [Variant(0, 1, ""), Variant(1, 2, "")], [Variant(0, 4, "")], Relation.IS_CONTAINED),
    ("CC", [Variant(0, 0, "C"), Variant(1, 2, "")], [Variant(0, 0, "C"), Variant(1, 1, "T")], Relation.DISJOINT),
    ("CTCCG", [Variant(4, 4, "TCC")], [Variant(0, 0, "GT")], Relation.OVERLAP),
    ("CCCCC", [Variant(4, 5, "")], [Variant(0, 1, "")], Relation.EQUIVALENT),
    ('AACCCCTTGTTC', [Variant(2, 2, "C")], [Variant(6, 7, "C")], Relation.IS_CONTAINED),
    ('AACATTTC', [Variant(1, 1, "AC")], [Variant(4, 5, "C")], Relation.OVERLAP),
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
