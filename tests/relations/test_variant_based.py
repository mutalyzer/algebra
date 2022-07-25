import pytest
from algebra import Relation, Variant, compare


@pytest.mark.parametrize("reference, lhs, rhs, expected", [
])
def test_compare(reference, lhs, rhs, expected):
    assert compare(reference, lhs, rhs) == expected
