import pytest

from algebra.extractor import extract
from algebra.variants import Parser, patch

variants = [
    ("CATATAGT CATAGAT", "5_7delinsGA"),
    ("CTAA TTA", "1_3delinsTT"),  # - with equals inside
    ("CTAACG TTACC", "1_6delinsTTACC"),  # - with equals inside
    ("CTTTG CTATTTT", "3_5delinsATTTT"),  # - with consecutive equals inside
    ("GCCTT GCAGCCCAT", "3_4delinsAGCCCA"),  # - with consecutive equals inside
    ("AGGTA AAGAAGGGGA", "2_4delinsAGAAGGGG"),  # - with consecutive equals inside
    ("TTGTA TTTGTGTT", "3_5delinsTGTGTT"),  # - with consecutive equals inside
    (
        "ACTAA ACGCCTATTAAATAAA",
        "3delinsGCCTATTAAATA",
    ),  # - with consecutive equals inside
    (
        "CAGGG AACTCAGGTAGGGTTAGAT",
        "1_5delinsAACTCAGGTAGGGTTAGAT",
    ),  # - with three consecutive equals inside
    (
        "CGACTGACGTTACCGAAGTTTTTTGTACAGTCGACTGACG CGACTGACATTACCGAAGTTTTTTTGTACAGGGTTCTGACG",
        "[9G>A;19_24T[7];31_34delinsGGTT]",
    ),  # - complex component
    (
        "CGACTGACGTTACCGAAGTTTTTTGTACAGTCGACTGACGTTCGTCCATGATACAGAGTATGCGCAATTCC CGACTGACATTACCGAAGTTTTTTTGTACAGGGTTCTGACGATCGTCCATGGCACGGGTATGCGCGCAATTGC",
        "[9G>A;19_24T[7];31_34delinsGGTT;41T>A;51_57delinsGCACGG;62_65GC[3];70C>G]",
    ),  # - complex component
    ("GTGCCCTAAGGGAT GAGCCTTAGGGCT", "[2T>A;6_8delinsTT;13A>C]"),
    ("TCATCAT TCAT", "2_7CAT[1]"),
    ("CATCAT CAT", "1_6CAT[1]"),
    ("CATCAT TCAT", "1_2del"),
]


@pytest.mark.parametrize("sequences, expected", variants)
def test_observed(sequences, expected):
    reference, observed = sequences.split(" ")
    new_variants = extract(reference, observed)

    assert patch(reference, Parser(new_variants).hgvs()) == observed


@pytest.mark.parametrize("sequences, expected", variants)
def test_variants(sequences, expected):
    reference, observed = sequences.split(" ")
    assert extract(reference, observed) == expected
