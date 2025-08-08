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
