import pytest
from algebra.lcs.all_lcs import _Node, edit


@pytest.mark.parametrize("reference, observed, expected_distance, expected_lcs_nodes", [
    ("", "", 0, []),
    ("AA", "ACA", 1, [[_Node(1, 1, 1)], [_Node(2, 3, 1)]]),
    ("ACA", "AA", 1, [[_Node(1, 1, 1)], [_Node(3, 2, 1)]]),
])
def test_edit(reference, observed, expected_distance, expected_lcs_nodes):
    distance, lcs_nodes = edit(reference, observed)
    assert distance == expected_distance
    assert lcs_nodes == expected_lcs_nodes
