import pytest
from algebra.lcs import LCSgraph
from algebra.variants import Variant


@pytest.mark.parametrize("args, expected", [
    (("GTGTGTTTTTTTAACAGGGA", [Variant(8, 9, "")]), [Variant(5, 12, "TTTTTT")]),
])
def test_canonical(args, expected):
    graph = LCSgraph.from_variants(*args)
    assert graph.canonical() == expected


@pytest.mark.parametrize("args, expected", [
    (("GTGTGTTTTTTTAACAGGGA", [Variant(8, 9, "")]), [Variant(5, 12, "TTTTTT")]),
])
def test_local_supremal(args, expected):
    graph = LCSgraph.from_variants(*args)
    assert graph.local_supremal() == expected


@pytest.mark.parametrize("args, expected", [
    (("GTGTGTTTTTTTAACAGGGA", [Variant(8, 9, "")]), Variant(5, 12, "TTTTTT")),
])
def test_supremal(args, expected):
    graph = LCSgraph.from_variants(*args)
    assert graph.supremal() == expected
