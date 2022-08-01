import pytest
from algebra import Relation
from algebra.relations.sequence_based import (are_disjoint, are_equivalent, compare,
                                              contains, have_overlap, is_contained)


TESTS = [
    ("A", "B", "B", Relation.EQUIVALENT),
    ("AAA", "AAB", "AAB", Relation.EQUIVALENT),
    ("AAA", "AA", "AA", Relation.EQUIVALENT),
    ("AAA", "CAAA", "CAAA", Relation.EQUIVALENT),
    ("AAA", "", "", Relation.EQUIVALENT),
    ("AAA", "ABB", "ABB", Relation.EQUIVALENT),
    ("AA", "AB", "BB", Relation.IS_CONTAINED),
    ("AAA", "AAB", "ABB", Relation.IS_CONTAINED),
    ("", "A", "AA", Relation.IS_CONTAINED),
    ("", "AB", "ABAB", Relation.IS_CONTAINED),
    ("ATATA", "ATATAB", "ATBTAB", Relation.IS_CONTAINED),
    ("", "BB", "BAB", Relation.IS_CONTAINED),
    ("CATATATC", "CATATTATC", "CATATATATC", Relation.IS_CONTAINED),
    ("AA", "BB", "AB", Relation.CONTAINS),
    ("AAA", "ABB", "AAB", Relation.CONTAINS),
    ("", "AA", "A", Relation.CONTAINS),
    ("", "ABAB", "AB", Relation.CONTAINS),
    ("ATATA", "ATBTAB", "ATATAB", Relation.CONTAINS),
    ("", "BAB", "BB", Relation.CONTAINS),
    ("CATATATC", "CATATATATC", "CATATTATC", Relation.CONTAINS),
    ("A", "B", "C", Relation.OVERLAP),
    ("AAA", "ABC", "ABD", Relation.OVERLAP),
    ("AAA", "BBA", "ABB", Relation.OVERLAP),
    ("", "BC", "CAB", Relation.OVERLAP),
    ("ATA", "BTA", "ATB", Relation.DISJOINT),
    ("AAA", "BAA", "AAA", Relation.DISJOINT),
    ("AAA", "AAA", "AAB", Relation.DISJOINT),
    ("AAA", "BAAA", "AAAB", Relation.DISJOINT),
    ("AAA", "AAAB", "BAAA", Relation.DISJOINT),
    ("", "A", "B", Relation.DISJOINT),
    ("T", "GG", "GGTA", Relation.OVERLAP),
    ("TC", "GTC", "GAA", Relation.IS_CONTAINED),
    ("T", "GC", "CT", Relation.CONTAINS),
    ("CT", "TG", "GC", Relation.DISJOINT),
    ("A", "ABD", "ABC", Relation.OVERLAP),
    ("A", "AB", "AC", Relation.DISJOINT),
    ("AAA", "BAAA", "AAAB", Relation.DISJOINT),
    ("A", "BAC", "BAD", Relation.OVERLAP),
    ("AA", "BAAC", "BAAD", Relation.OVERLAP),
    ("AAA", "BAAAC", "BAAAD", Relation.OVERLAP),
    ("TGTA", "CTGCT", "TAGGAACG", Relation.DISJOINT),
    ("CT", "GT", "AT", Relation.OVERLAP),
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
    relation = compare(reference, lhs, rhs)
    assert relation == expected
