import pytest

from algebra.extractor.extractor_dfa import extract


@pytest.mark.parametrize(
    "reference, observed, variants",
    [
        ("TAA", "T", ["2_3del"]),
        ("TAACCAG", "GATCC", ["1_3delinsGAT", "6_7del"]),
        ("CGTGA", "CCGGATATT", ["1delinsCC", "3del", "5_6insTATT"]),
        ("GAAG", "GGAAGCACG", ["1delinsGG", "4delinsGCACG"]),
        ("GAAGC", "GGAAGCACGC", ["1delinsGG", "4_5delinsGCACGC"]),
        ("CTGAAT", "TCAAATTG", ["1_3delinsTCA", "6_7insTG"]),
        ("AGATAGCCTAACGT", "AGCCT", ["1_6delinsAG", "9_14delinsT"]),
        ("GAGTTA", "AGGTATG", ["1del", "4T>G", "6_7insTG"]),

    ],
)
def test_extractor_mdfa(reference, observed, variants):
    assert [e.to_hgvs(reference) for e in extract(reference, observed)] == variants
