import pytest
from algebra import Variant
from algebra.extractor.local_suprema import local_suprema
from algebra.lcs import lcs_graph


@pytest.mark.parametrize(
    "reference, observed, expected_suprema, hgvs",
    [
        (
            "",
            "",
            [],
            "",
        ),
        (
            "",
            "C",
            [Variant(0, 0, "C")],
            "0_1insC",
        ),
        (
            "C",
            "",
            [Variant(0, 1, "")],
            "1del",
        ),
        (
            "C",
            "C",
            [],
            "",
        ),
        (
            "C",
            "G",
            [Variant(0, 1, "G")],
            "1C>G",
        ),
        (
            "CATC",
            "GATG",
            [Variant(0, 1, "G"),
             Variant(3, 4, "G")],
            "1C>G;4C>G",
        ),
        (
            "AGAATTGCTTGAA",
            "AGGGTTAAA",
            [Variant(1, 8, "GGG"), Variant(10, 13, "AAA")],
            "2_8delinsGGG;11_13delinsAAA",
        ),
        (
            "CTCTAGAGACTTTATTTTCCAC",
            "GTCTCAGACTTTCTTTATCCCC",
            [Variant(0, 9, "GTCTCAGA"),
             Variant(13, 14, "C"),
             Variant(17, 17, "A"),
             Variant(18, 22, "CCCC")
             ],
            "1_9delinsGTCTCAGA;14A>C;17_18insA;19_22delinsCCCC",
        ),
        (
            "ATATACCTTTTA",
            "CTATAGCCTTTTTC",
            [Variant(0, 1, "C"),
             Variant(5, 5, "G"),
             Variant(7, 12, "TTTTTC")
             ],
            "1A>C;5_6insG;8_12delinsTTTTTC",
        ),
        (
            "CAGGGGAAGTG",
            "GCAGGGGCCTA",
            [Variant(0, 0, "G"),
             Variant(2, 11, "GGGGCCTA")
             ],
            "0_1insG;3_11delinsGGGGCCTA",
        ),
        (
            "TCGTGGT",
            "CTAACAT",
            [Variant(0, 7, "CTAACAT")],
            "1_7delinsCTAACAT",
        ),
    ],
)
def test_local_suprema(reference, observed, expected_suprema, hgvs):
    suprema = list(local_suprema(reference, observed, lcs_graph(reference, observed)))
    assert suprema == expected_suprema
    assert ";".join([supremum.to_hgvs(reference) for supremum in suprema]) == hgvs
