import pytest
from algebra import Variant
from algebra.lcs.all_lcs import (_Node, bfs_traversal, build_graph,
                                 dfs_traversal, edit, lcs_graph)


@pytest.mark.parametrize("reference, observed, expected_distance, expected_lcs_nodes", [
    ("", "", 0, []),
    ("AA", "ACA", 1, [[_Node(0, 0, 1)], [_Node(1, 2, 1)]]),
    ("ACA", "AA", 1, [[_Node(0, 0, 1)], [_Node(2, 1, 1)]]),
    ("CTCGGCATTA", "GGCTGGCTGT", 6, [
        [_Node(2, 2, 1), _Node(3, 1, 1)],
        [_Node(0, 2, 2)],
        [_Node(4, 4, 1), _Node(3, 0, 3), _Node(3, 5, 1)],
        [],
        [_Node(3, 4, 3)],
        [_Node(7, 7, 1), _Node(8, 7, 1)],
        [_Node(8, 9, 1)]
    ]),
    ("CATATATCG", "CTTATAGCAT", 7, [
        [_Node(0, 0, 1)],
        [_Node(2, 1, 1), _Node(4, 1, 1), _Node(1, 5, 1)],
        [],
        [_Node(1, 3, 3)],
        [_Node(2, 2, 4), _Node(4, 2, 3)],
        [_Node(7, 7, 1), _Node(8, 6, 1), _Node(5, 8, 2)],
    ]),
    ("TTT", "TTTTAT", 3, [
        [_Node(0, 3, 1)],
        [_Node(0, 2, 2)],
        [_Node(0, 0, 3), _Node(0, 1, 3), _Node(2, 5, 1)],
    ]),
])
def test_edit(reference, observed, expected_distance, expected_lcs_nodes):
    distance, lcs_nodes = edit(reference, observed)
    assert distance == expected_distance
    assert distance == len(reference) - len(lcs_nodes) + len(observed) - len(lcs_nodes)
    assert lcs_nodes == expected_lcs_nodes


@pytest.mark.parametrize("reference, observed, max_distance, expected_distance", [
    ("CTCGGCATTA", "GGCTGGCTGT", 6, 6),
])
def test_edit_max_distance(reference, observed, max_distance, expected_distance):
    distance, _ = edit(reference, observed, max_distance=max_distance)
    assert distance == expected_distance


@pytest.mark.parametrize("reference, observed, max_distance, exception, message", [
    ("CTCGGCATTA", "GGCTGGCTGT", 5, ValueError, "maximum distance exceeded"),
])
def test_edit_max_distance_fail(reference, observed, max_distance, exception, message):
    with pytest.raises(exception) as exc:
        distance, _ = edit(reference, observed, max_distance=max_distance)
    assert str(exc.value) == message


@pytest.mark.parametrize("reference, observed, expected_edges", [
    ("", "", []),
    ("A", "A", []),
    ("TTAATTGACA", "CTACTGAGTT", {
        Variant(8, 10, "GTT"), Variant(10, 10, "GTT"), Variant(6, 10, ""),
        Variant(7, 9, ""),  Variant(3, 4, "G"), Variant(4, 4, "G"),
        Variant(5, 6, ""), Variant(3, 4, "C"), Variant(4, 4, "C"),
        Variant(3, 5, "C"), Variant(4, 5, "C"), Variant(2, 2, "G"),
        Variant(3, 3, "CTG"), Variant(2, 3, "G"),  Variant(1, 2, ""),
        Variant(1, 3, ""), Variant(2, 3, ""), Variant(1, 1, "AC"),
        Variant(0, 0, "C"), Variant(0, 1, "C")
    }),
    ("TTT", "TATTTT", {
        Variant(3, 3, "TT"), Variant(3, 3, "T"), Variant(2, 2, "T"),
        Variant(2, 2, "TT"), Variant(1, 1, "A"),
        Variant(1, 1, "AT"), Variant(1, 1, "ATT"), Variant(1, 1, "T"),
        Variant(0, 0, "TA"), Variant(0, 0, "TAT")
    }),
    ("TCTCTATCGTA", "TCTA", {
        Variant(6, 11, ""), Variant(3, 5, ""), Variant(3, 10, ""),
        Variant(5, 10, ""), Variant(7, 10, ""), Variant(2, 4, ""),
        Variant(2, 6, ""), Variant(4, 6, ""), Variant(2, 9, ""),
        Variant(4, 9, ""), Variant(8, 9, ""), Variant(1, 3, ""),
        Variant(1, 7, ""), Variant(3, 7, ""), Variant(5, 7, ""),
        Variant(0, 2, ""), Variant(0, 4, ""), Variant(0, 6, ""),
    }),
    ("AAAATA", "GAAAAGAAA", {
        Variant(6, 6, "AA"), Variant(6, 6, "A"), Variant(4, 5, "G"),
        Variant(4, 5, "GA"), Variant(4, 5, ""), Variant(4, 5, "GAA"),
        Variant(4, 5, "A"), Variant(3, 3, "AG"),
        Variant(3, 3, "G"), Variant(3, 3, "AGA"), Variant(3, 3, "GA"),
        Variant(2, 2, "A"), Variant(2, 2, "AAG"), Variant(2, 2, "AG"),
        Variant(2, 2, "G"), Variant(1, 1, "A"), Variant(1, 1, "AA"),
        Variant(1, 1, "A"), Variant(0, 0, "G"), Variant(0, 0, "GA"),
        Variant(0, 0, "GAA")
    }),
    ("CATATATCG", "CTTATAGCAT", {
        Variant(1, 1, "TT"), Variant(1 ,2, ""), Variant(1, 2, "T"),
        Variant(4, 5, "GC"), Variant(3, 3, "T"), Variant(3, 4, ""),
        Variant(6, 7, "G"), Variant(6, 8, ""), Variant(6, 6, "GCA"),
        Variant(5, 5, "AGC"), Variant(7, 9, ""), Variant(7, 7, "AG"),
        Variant(7, 8, "A"), Variant(8, 9, "AT"), Variant(9, 9, "CAT"),
    }),
    ("TTT", "TTTTAT", {
        # TODO
    }),
])
def test_build_graph(reference, observed, expected_edges):
    _, lcs_nodes = edit(reference, observed)
    _, edges = build_graph(reference, observed, lcs_nodes)
    assert edges == expected_edges


@pytest.mark.parametrize("reference, observed, edges", [
    ("", "", []),
    ("TTGGTT", "AAGTTAG", [
        [Variant(0, 0, "AAG")],
        [Variant(0, 2, "AA")],
        [Variant(0, 3, "AA")],
        [Variant(2, 2, "A")],
        [Variant(2, 3, "A")],
        [Variant(3, 4, "")],
        [Variant(6, 6, "AG")],
        [Variant(3, 6, "")],
        [Variant(4, 6, "")],
     ]),
])
def test_bfs_traversal(reference, observed, edges):
    graph = lcs_graph(reference, observed)
    assert [edge for *_, edge in bfs_traversal(graph)] == edges


@pytest.mark.parametrize("reference, observed, edges", [
    ("", "", []),
    ("TTGGTT", "AAGTTAG", [
        [Variant(0, 0, "AAG")],
        [Variant(0, 0, "AA"), Variant(0, 1, ""), Variant(1, 2, "")],
        [Variant(0, 0, "A"), Variant(0, 1, ""), Variant(1, 1, "A"), Variant(1, 2, "")],
        [Variant(0, 0, "A"), Variant(0, 1, ""), Variant(1, 2, ""), Variant(2, 2, "A")],
        [Variant(0, 1, ""), Variant(1, 1, "AA"), Variant(1, 2, "")],
        [Variant(0, 1, ""), Variant(1, 1, "A"), Variant(1, 2, ""), Variant(2, 2, "A")],
        [Variant(0, 1, ""), Variant(1, 2, ""), Variant(2, 2, "AA")],
        [Variant(0, 0, "AA"), Variant(0, 1, ""), Variant(1, 2, ""), Variant(2, 3, "")],
        [Variant(0, 0, "A"), Variant(0, 1, ""), Variant(1, 1, "A"), Variant(1, 2, ""), Variant(2, 3, "")],
        [Variant(0, 0, "A"), Variant(0, 1, ""), Variant(1, 2, ""), Variant(2, 2, "A"), Variant(2, 3, "")],
        [Variant(0, 0, "A"), Variant(0, 1, ""), Variant(1, 2, ""), Variant(2, 3, ""), Variant(3, 3, "A")],
        [Variant(0, 1, ""), Variant(1, 1, "AA"), Variant(1, 2, ""), Variant(2, 3, "")],
        [Variant(0, 1, ""), Variant(1, 1, "A"), Variant(1, 2, ""), Variant(2, 2, "A"), Variant(2, 3, "")],
        [Variant(0, 1, ""), Variant(1, 1, "A"), Variant(1, 2, ""), Variant(2, 3, ""), Variant(3, 3, "A")],
        [Variant(0, 1, ""), Variant(1, 2, ""), Variant(2, 2, "AA"), Variant(2, 3, "")],
        [Variant(0, 1, ""), Variant(1, 2, ""), Variant(2, 2, "A"), Variant(2, 3, ""), Variant(3, 3, "A")],
        [Variant(0, 1, ""), Variant(1, 2, ""), Variant(2, 3, ""), Variant(3, 3, "AA")],
        [Variant(2, 2, "A")],
        [Variant(2, 2, "A"), Variant(2, 3, "")],
        [Variant(2, 3, ""), Variant(3, 3, "A")],
        [Variant(3, 4, "")],
        [Variant(6, 6, "AG")],
        [Variant(3, 4, ""), Variant(4, 5, ""), Variant(5, 6, "")],
        [Variant(4, 5, ""), Variant(5, 6, "")],
     ]),
])
def test_bfs_traversal_atomics(reference, observed, edges):
    graph = lcs_graph(reference, observed)
    assert [edge for *_, edge in bfs_traversal(graph, atomics=True)] == edges


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
def test_dfs_traversal(reference, observed, expected_variant):
    graph = lcs_graph(reference, observed)
    assert list(dfs_traversal(graph)) == expected_variant


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
def test_dfs_traversal_atomics(reference, observed, expected_variant):
    graph = lcs_graph(reference, observed)
    assert list(dfs_traversal(graph, atomics=True)) == expected_variant


@pytest.mark.parametrize("duplicates, unique", [
    ([_Node(0, 0), _Node(0, 0), _Node(0, 0)], {_Node(0, 0)}),
])
def test_node_hash(duplicates, unique):
    assert set(duplicates) == unique
