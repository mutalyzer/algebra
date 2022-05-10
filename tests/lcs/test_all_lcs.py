import pytest
from algebra.lcs import edit, lcs_graph, maximal_variant, traversal
from algebra.lcs.all_lcs import _Node
from algebra.variants import Variant


@pytest.mark.parametrize("reference, observed, expected_distance, expected_lcs_nodes", [
    ("", "", 0, []),
    ("AA", "ACA", 1, [[_Node(1, 1, 1)], [_Node(2, 3, 1)]]),
    ("ACA", "AA", 1, [[_Node(1, 1, 1)], [_Node(3, 2, 1)]]),
    ("CTCGGCATTA", "GGCTGGCTGT", 6, [
        [_Node(3, 3, 1), _Node(4, 2, 1)],
        [_Node(1, 3, 2)],
        [_Node(5, 5, 1), _Node(4, 1, 3), _Node(4, 6, 1)],
        [],
        [_Node(4, 5, 3)],
        [_Node(8, 8, 1), _Node(9, 8, 1)],
        [_Node(9, 10, 1)]
    ]),
])
def test_edit(reference, observed, expected_distance, expected_lcs_nodes):
    distance, lcs_nodes = edit(reference, observed)
    for level in lcs_nodes:
        print(level)
    assert distance == expected_distance
    assert lcs_nodes == expected_lcs_nodes


@pytest.mark.parametrize("reference, observed, expected_edges", [
    ("", "", []),
    ("TTAATTGACA", "CTACTGAGTT", [
        Variant(8, 10, "GTT"), Variant(10, 10, "GTT"), Variant(6, 10),
        Variant(7, 9), Variant(5, 6), Variant(3, 4, "G"), Variant(4, 4, "G"),
        Variant(3, 4, "C"), Variant(4, 4, "C"), Variant(2, 2, "G"),
        Variant(3, 3, "CTG"), Variant(2, 3, "G"), Variant(3, 5, "C"),
        Variant(4, 5, "C"), Variant(1, 2), Variant(1, 3), Variant(2, 3),
        Variant(1, 1, "AC"), Variant(0, 0, "C"), Variant(0, 1, "C")
    ]),
    ("TTT", "TATTTT", [
        Variant(3, 3, "TT"), Variant(3, 3, "T"), Variant(2, 2, "T"),
        Variant(2, 2, "TT"), Variant(2, 2, "T"), Variant(1, 1, "A"),
        Variant(1, 1, "AT"), Variant(1, 1, "ATT"), Variant(1, 1, "T"),
        Variant(0, 0, "TA"), Variant(0, 0, "TAT")
    ]),
    ("TCTCTATCGTA", "TCTA", [
        Variant(6, 11), Variant(3, 5), Variant(3, 10), Variant(7, 10),
        Variant(5, 10), Variant(4, 6), Variant(2, 6), Variant(2, 4),
        Variant(8, 9), Variant(2, 9), Variant(4, 9), Variant(1, 7),
        Variant(3, 7), Variant(5, 7), Variant(1, 3), Variant(0, 4),
        Variant(0, 6), Variant(0, 2)
    ]),
    ("AAAATA", "GAAAAGAAA", [
        Variant(6, 6, "AA"), Variant(6, 6, "A"), Variant(4, 5, "G"),
        Variant(4, 5, "GA"), Variant(4, 5), Variant(4, 5, "GAA"),
        Variant(4, 5, "A"), Variant(4, 5), Variant(3, 3, "G"),
        Variant(3, 3, "AG"), Variant(3, 3, "GA"), Variant(3, 3, "AGA"),
        Variant(2, 2, "A"), Variant(2, 2, "G"), Variant(2, 2, "AG"),
        Variant(2, 2, "AAG"), Variant(1, 1, "A"), Variant(1, 1, "AA"),
        Variant(1, 1, "A"), Variant(0, 0, "GAA"), Variant(0, 0, "GA"),
        Variant(0, 0, "G")
    ]),
])
def test_lcs_graph(reference, observed, expected_edges):
    _, lcs_nodes = edit(reference, observed)
    _, edges = lcs_graph(reference, observed, lcs_nodes)
    assert edges == expected_edges


@pytest.mark.parametrize("reference, observed, expected_variant", [
    ("", "", Variant(0, 0)),
    ("A", "C", Variant(0, 1, "C")),
    ("", "A", Variant(0, 0, "A")),
    ("A", "", Variant(0, 1)),
    ("ACCCA", "ACCA", Variant(1, 4, "CC")),
    ("AACCT", "AACGT", Variant(2, 4, "CG")),
])
def test_lcs_graph_max_variant(reference, observed, expected_variant):
    _, lcs_nodes = edit(reference, observed)
    _, edges = lcs_graph(reference, observed, lcs_nodes)
    assert maximal_variant(reference, observed, edges) == expected_variant


@pytest.mark.parametrize("reference, observed, expected_variant", [
    ("", "", [[]]),
    ("TTT", "TATTTT", [
        [Variant(1, 1, "A"), Variant(3, 3, "TT")],
        [Variant(1, 1, "A"), Variant(2, 2, "T"), Variant(3, 3, "T")],
        [Variant(1, 1, "A"), Variant(2, 2, "TT")],
        [Variant(1, 1, "AT"), Variant(3, 3, "T")],
        [Variant(1, 1, "AT"), Variant(2, 2, "T")],
        [Variant(1, 1, "ATT")],
        [Variant(0, 0, "TA"), Variant(3, 3, "T")],
        [Variant(0, 0, "TA"), Variant(2, 2, "T")],
        [Variant(0, 0, "TA"), Variant(1, 1, "T")],
        [Variant(0, 0, "TAT")]
    ]),
])
def test_traversal(reference, observed, expected_variant):
    _, lcs_nodes = edit(reference, observed)
    root, _ = lcs_graph(reference, observed, lcs_nodes)
    assert list(traversal(root)) == expected_variant


@pytest.mark.parametrize("reference, observed, expected_variant", [
    ("", "", [[]]),
    ("TTT", "TATTTT", [
        [Variant(1, 1, "A"), Variant(3, 3, "TT")],
        [Variant(1, 1, "A"), Variant(2, 2, "T"), Variant(3, 3, "T")],
        [Variant(1, 1, "A"), Variant(2, 2, "TT")],
        [Variant(1, 1, "AT"), Variant(3, 3, "T")],
        [Variant(1, 1, "AT"), Variant(2, 2, "T")],
        [Variant(1, 1, "ATT")],
        [Variant(0, 0, "TA"), Variant(3, 3, "T")],
        [Variant(0, 0, "TA"), Variant(2, 2, "T")],
        [Variant(0, 0, "TA"), Variant(1, 1, "T")],
        [Variant(0, 0, "TAT")]
    ]),
])
def test_traversal_atomics(reference, observed, expected_variant):
    _, lcs_nodes = edit(reference, observed)
    root, _ = lcs_graph(reference, observed, lcs_nodes)
    assert list(traversal(root, atomics=True)) == expected_variant
