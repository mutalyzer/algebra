import pytest
from algebra.lcs import LCSgraph
from algebra.variants import Variant


@pytest.mark.parametrize("reference, variants, supremal, local_supremal, canonical", [
    ("GTGTGTTTTTTTAACAGGGA", [Variant(8, 9, "")], Variant(5, 12, "TTTTTT"), [Variant(5, 12, "TTTTTT")], [Variant(5, 12, "TTTTTT")]),
    ("GTGTGTTTTTTTAACAGGGA", [Variant(4, 5, ""), Variant(12, 13, "")], Variant(4, 14, "TTTTTTTA"), [Variant(4, 5, ""), Variant(12, 14, "A")], [Variant(4, 5, ""), Variant(12, 14, "A")]),
    ("GGGTCTTCGACTTTCCACGAAAATCGC", [Variant(11, 12, "")], Variant(11, 14, "TT"), [Variant(11, 14, "TT")], [Variant(11, 14, "TT")]),
])
def test_LCSgraph_from_variants(reference, variants, supremal, local_supremal, canonical):
    graph = LCSgraph.from_variants(reference, variants)
    assert all((graph.supremal() == supremal,
        graph.local_supremal() == local_supremal,
        graph.canonical() == canonical))


def test_LCSgraph_edges():
    graph = LCSgraph.from_variants("ACCTGACT", [Variant(1, 2, "T"), Variant(4, 5, "T"), Variant(7, 8, "TT")])
    assert list(graph.edges()) == [
        {"head": (5, 5, 3), "tail": (7, 8, 1), "variant": Variant(7, 7, "T"), "count": 2},
        {"head": (3, 3, 1), "tail": (5, 5, 3), "variant": Variant(4, 5, "T"), "count": 1},
        {"head": (3, 4, 1), "tail": (5, 5, 3), "variant": Variant(4, 5, ""), "count": 1},
        {"head": (1, 2, 1), "tail": (3, 4, 1), "variant": Variant(2, 3, "T"), "count": 1},
        {"head": (1, 2, 1), "tail": (3, 3, 1), "variant": Variant(2, 3, ""), "count": 1},
        {"head": (2, 2, 1), "tail": (3, 3, 1), "variant": None, "count": 0},
        {"head": (2, 2, 1), "tail": (3, 4, 1), "variant": Variant(3, 3, "T"), "count": 1},
        {"head": (1, 1, 0), "tail": (2, 2, 1), "variant": Variant(1, 2, "T"), "count": 1},
        {"head": (1, 1, 0), "tail": (1, 2, 1), "variant": Variant(1, 1, "T"), "count": 1}
    ]
