import pytest
from test_variant_based import TESTS
from algebra import LCSgraph, Relation
from algebra.relations.graph_based import (are_disjoint, are_equivalent, compare, contains,
                                           have_overlap, is_contained)


@pytest.mark.parametrize("reference, lhs, rhs, expected", TESTS)
def test_are_equivalent(reference, lhs, rhs, expected):
    assert are_equivalent(reference, LCSgraph.from_variant(reference, lhs), LCSgraph.from_variant(reference, rhs)) == (expected == Relation.EQUIVALENT)


@pytest.mark.parametrize("reference, lhs, rhs, expected", TESTS)
def test_contains(reference, lhs, rhs, expected):
    assert contains(reference, LCSgraph.from_variant(reference, lhs), LCSgraph.from_variant(reference, rhs)) == (expected == Relation.CONTAINS)


@pytest.mark.parametrize("reference, lhs, rhs, expected", TESTS)
def test_is_contained(reference, lhs, rhs, expected):
    assert is_contained(reference, LCSgraph.from_variant(reference, lhs), LCSgraph.from_variant(reference, rhs)) == (expected == Relation.IS_CONTAINED)


@pytest.mark.parametrize("reference, lhs, rhs, expected", TESTS)
def test_have_overlap(reference, lhs, rhs, expected):
    assert have_overlap(reference, LCSgraph.from_variant(reference, lhs), LCSgraph.from_variant(reference, rhs)) == (expected == Relation.OVERLAP)


@pytest.mark.parametrize("reference, lhs, rhs, expected", TESTS)
def test_are_disjoint(reference, lhs, rhs, expected):
    assert are_disjoint(reference, LCSgraph.from_variant(reference, lhs), LCSgraph.from_variant(reference, rhs)) == (expected == Relation.DISJOINT)


@pytest.mark.parametrize("reference, lhs, rhs, expected", TESTS)
def test_compare(reference, lhs, rhs, expected):
    assert compare(reference, LCSgraph.from_variant(reference, lhs), LCSgraph.from_variant(reference, rhs)) == expected
