import pytest
from algebra import Variant
from algebra.lcs.all_lcs import edit
from algebra.lcs.all_lcs_dfa import _Node
from algebra.lcs.all_lcs_dfa import lcs_graph_dfa as lcs_graph
from algebra.lcs.all_lcs_dfa import traversal


@pytest.mark.parametrize("reference, observed, expected_edges", [
    ("", "", []),
    ("A", "A", []),
    ("TTAATTGACA", "CTACTGAGTT", [
        Variant(8, 10, "GTT"), Variant(10, 10, "GTT"), Variant(6, 10, ""),
        Variant(7, 9, ""), Variant(5, 6, ""), Variant(3, 4, "G"), Variant(4, 4, "G"),
        Variant(3, 4, "C"), Variant(4, 4, "C"), Variant(2, 2, "G"),
        Variant(3, 3, "CTG"), Variant(2, 3, "G"), Variant(3, 5, "C"),
        Variant(4, 5, "C"), Variant(1, 2, ""), Variant(1, 3, ""), Variant(2, 3, ""),
        Variant(1, 1, "AC"), Variant(0, 0, "C"), Variant(0, 1, "C")
    ]),
    ("TTT", "TATTTT", [
        Variant(3, 3, "TT"), Variant(3, 3, "T"), Variant(2, 2, "T"),
        Variant(2, 2, "TT"), Variant(2, 2, "T"), Variant(1, 1, "A"),
        Variant(1, 1, "AT"), Variant(1, 1, "ATT"), Variant(1, 1, "T"),
        Variant(0, 0, "TA"), Variant(0, 0, "TAT")
    ]),
    ("TCTCTATCGTA", "TCTA", [
        Variant(6, 11, ""), Variant(3, 5, ""), Variant(3, 10, ""), Variant(7, 10, ""),
        Variant(5, 10, ""), Variant(4, 6, ""), Variant(2, 6, ""), Variant(2, 4, ""),
        Variant(8, 9, ""), Variant(2, 9, ""), Variant(4, 9, ""), Variant(1, 7, ""),
        Variant(3, 7, ""), Variant(5, 7, ""), Variant(1, 3, ""), Variant(0, 4, ""),
        Variant(0, 6, ""), Variant(0, 2, "")
    ]),
    ("AAAATA", "GAAAAGAAA", [
        Variant(6, 6, "AA"), Variant(6, 6, "A"), Variant(4, 5, "G"),
        Variant(4, 5, "GA"), Variant(4, 5, ""), Variant(4, 5, "GAA"),
        Variant(4, 5, "A"), Variant(4, 5, ""), Variant(3, 3, "G"),
        Variant(3, 3, "AG"), Variant(3, 3, "GA"), Variant(3, 3, "AGA"),
        Variant(2, 2, "A"), Variant(2, 2, "G"), Variant(2, 2, "AG"),
        Variant(2, 2, "AAG"), Variant(1, 1, "A"), Variant(1, 1, "AA"),
        Variant(1, 1, "A"), Variant(0, 0, "GAA"), Variant(0, 0, "GA"),
        Variant(0, 0, "G")
    ]),
])
def test_lcs_graph(reference, observed, expected_edges):
    print(reference, observed)
    _, lcs_nodes = edit(reference, observed)
    _, edges = lcs_graph(reference, observed, lcs_nodes)
    print(lcs_nodes)
    print(edges)
    assert sorted(edges, key=lambda e: (e.start, e.end, e.sequence)) == sorted(expected_edges, key=lambda e: (e.start, e.end, e.sequence))


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
    print(lcs_nodes)
    print(root)
    print(root.edges)
    print(list(traversal(root)))
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


@pytest.mark.parametrize("duplicates, unique", [
    ([_Node(0, 0), _Node(0, 0), _Node(0, 0)], [_Node(0, 0)]),
])
def test_node_hash(duplicates, unique):
    assert set(duplicates) == set(unique)


@pytest.mark.parametrize("node, string", [
    (_Node(0, 0), "(0, 0, 0)"),
])
def test_node_string(node, string):
    assert str(node) == string
