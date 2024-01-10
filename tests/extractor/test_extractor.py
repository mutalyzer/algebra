import pytest
from algebra import LCSgraph, Variant
from algebra.extractor.extractor import (diagonal, extract, extract_sequence,
                                         extract_supremal, to_hgvs)


@pytest.mark.parametrize("reference, observed, variants", [
    ("CATATATCG", "CTTATAGCAT", [Variant(1, 2, "T"), Variant(6, 7, "G"), Variant(8, 9, "AT")]),
])
def test_diagonal(reference, observed, variants):
    assert diagonal(reference, observed, LCSgraph.from_sequence(reference, observed)) == variants


@pytest.mark.parametrize("reference, observed, variants, hgvs", [
    ("AAATAATATAATAATTTAT", "AAATAATAATTTAT", [Variant(2, 13, "ATAATA")], "3_12ATAAT[1]"),
    ("ACATATATATTATATATAAAA", "ACATATATATAAAA", [Variant(3, 17, "TATATAT")], "4_17TATATAT[1]"),
    ("TTTCTCTCTCTCTCTC", "TTTCTCTCTCTCCTCTCTCTC", [Variant(7, 12, "CTCTCCTCTC")], "8_12dup"),
    ("AAATCTGTATGGAG", "AAATCTGTTGTTGTATGGAG", [Variant(5, 8, "TGTTGTTGT")], "6_8TGT[3]"),
    ("AAAGCATTTTAAA", "AAAGCATGTGTTTTAAA", [Variant(6, 7, "TGTGT")], "7_8insGT[2]"),
    ("CAATAAATACAG", "CAATACAG", [Variant(1, 9, "AATA")], "2_9AATA[1]"),
    ("AGTGCTTTGTTTTGTTATAATTAAC", "AGTGCTTTGTTATAATTAAC", [Variant(5, 16, "TTTGTT")], "6_15TTTGT[1]"),
    ("ATTGCATTTCTTCAATACTAATTTCTAAAGCCTTT", "ATTGCATTTCTTCAATACTAATTTCTAAAGCCTTTTTCTTCAATACTAATTTCTAAAGCCTTT", [Variant(6, 35, "TTTCTTCAATACTAATTTCTAAAGCCTTTTTCTTCAATACTAATTTCTAAAGCCTTT")], "8_35dup"),
    ("TTCTTCAATACTAATTTCTAAAGCCTTTCCG", "TTCTTCAATACTAATTTCTAAAGCCTTTTTCTTCAATACTAATTTCTAAAGCCTTTCCG", [Variant(0, 28, "TTCTTCAATACTAATTTCTAAAGCCTTTTTCTTCAATACTAATTTCTAAAGCCTTT")], "1_28dup"),
    ("AAGTCTCATGGCTATTTGCA", "AAGTCTCATGGCTATATGGCTATTTGCA", [Variant(7, 15, "ATGGCTATATGGCTAT")], "8_15dup"),
    ("AAAGGGAGAGAAGACT", "AAAGGGAGAGACT", [Variant(8, 14, "AGA")], "9_14AGA[1]"),
    ("AAAGGGAGAGACT", "AAAGGGAGAGAAGACT", [Variant(8, 11, "AGAAGA")], "9_11dup"),
    ("AAAGGGAGAGAAGAACT", "AAAGGGAGAGACT", [Variant(10, 15, "A")], "12_15del"),
    ("AAAGGGAGAGAAGACT", "AAAGGGAGAGCT", [Variant(10, 14, "")], "11_14del"),
    ("ATTTCCCACTGAAAAATAAATCCCACCGGGC", "ATTTCCACCGGGC", [Variant(4, 24, "CC")], "7_24del"),
    ("ATAT", "ATATATATA", [Variant(4, 4, "ATATA")], "4_5ins[AT[2];A]"),
    ("CAAAAAA", "TTTTTTG", [Variant(0, 7, "TTTTTTG")], "1_7inv"),
    ("AACTCAGGTAGGGTTAGAT", "CAGGG", [Variant(0, 19, "CAGGG")], "1_19delinsCAGGG"),
    ("CAGGG", "AACTCAGGTAGGGTTAGAT", [Variant(0, 5, "AACTCAGGTAGGGTTAGAT")], "1_5delinsAACTCAGGTAGGGTTAGAT"),
    ("GTGCCCTAAGGGAT", "GAGCCTTAGGGCT", [Variant(1, 2, "A"), Variant(3, 9, "CCTTA"), Variant(12, 13, "C")], "[2T>A;6_8delinsT[2];13A>C]"),
    ("CATCAT", "TCAT", [Variant(0, 2, "")], "1_2del"),
    ("ACTAA", "ACGCCTATTAAATAAA", [Variant(1, 5, "CGCCTATTAAATAAA")], "3delinsGCCTATTAAATA"),
    ("TTGTA", "TTTGTGTT", [Variant(0, 5, "TTTGTGTT")], "3_5delinsTGTGTT"),
    ("AGGTA", "AAGAAGGGGA", [Variant(0, 4, "AAGAAGGGG")], "2_4delinsAGAAGGGG"),
    ("GCCTT", "GCAGCCCAT", [Variant(0, 4, "GCAGCCCA")], "3_4delinsAGCCCA"),
    ("CTAACG", "TTACC", [Variant(0, 6, "TTACC")], "1_6delinsTTACC"),
    ("CTAA", "TTA", [Variant(0, 4, "TTA")], "1_3delinsT[2]"),
    ("CATATAGT", "CATAGAT", [Variant(1, 7, "ATAGA")], "5_7delinsGA"),
    ("CGC", "CATC", [Variant(1, 2, "AT")], "2delinsAT"),
    ("CGC", "CATATATC", [Variant(1, 2, "ATATAT")], "2delinsAT[3]"),
    ("CGC", "CATATAC", [Variant(1, 2, "ATATA")], "2delins[AT[2];A]"),
    ("CGGC", "CATATAC", [Variant(1, 3, "ATATA")], "2_3delins[AT[2];A]"),
    ("CAT", "CAT", [], "="),
    ("TA", "TAA", [Variant(1, 2, "AA")], "2dup"),
    ("CATCAT", "CATCAGGGGGGGT", [Variant(5, 5, "GGGGGGG")], "5_6insG[7]"),
    ("CATCA", "CATCAGGGGGGG", [Variant(5, 5, "GGGGGGG")], "5_6insG[7]"),
    ("CATCAT", "CATCATCAT", [Variant(0, 6, "CATCATCAT")], "1_6CAT[3]"),
    ("CATCAT", "CATCATCATCAT", [Variant(0, 6, "CATCATCATCAT")], "1_6CAT[4]"),
    ("AA", "AAA", [Variant(0, 2, "AAA")], "1_2A[3]"),
    ("AA", "AAAA", [Variant(0, 2, "AAAA")], "1_2A[4]"),
    ("CATCAT", "CAT", [Variant(0, 6, "CAT")], "1_6CAT[1]"),
    ("CATCAT", "", [Variant(0, 6, "")], "1_6del"),
    ("", "CATCAT", [Variant(0, 0, "CATCAT")], "0_1insCAT[2]"),
    ("CATCATC", "CATCATCATC", [Variant(0, 7, "CATCATCATC")], "1_6CAT[3]"),
    ("CATCA", "CATCATCAT", [Variant(5, 5, "TCAT")], "5_6insTCAT"),
    ("ATCAT", "CATCATCAT", [Variant(0, 0, "CATC")], "0_1insCATC"),
    ("TT", "TATT", [Variant(0, 1, "TAT")], "1_2insAT"),
    ("TATT", "TT", [Variant(0, 3, "T")], "2_3del"),
    ("CTTTG", "CTATTTT", [Variant(1, 5, "TATTTT")], "3_5delinsATTTT"),
    ("TTT", "T", [Variant(0, 3, "T")], "1_3T[1]"),
    ("TT", "T", [Variant(0, 2, "T")], "1_2T[1]"),
    ("T", "TTT", [Variant(0, 1, "TTT")], "1T[3]"),
    ("AAA", "AAAA", [Variant(0, 3, "AAAA")], "1_3A[4]"),
    ("AAA", "AAAAA", [Variant(0, 3, "AAAAA")], "1_3A[5]"),
    ("AAA", "AAAAAA", [Variant(0, 3, "AAAAAA")], "1_3A[6]"),
    ("AAA", "AAAAAAA", [Variant(0, 3, "AAAAAAA")], "1_3A[7]"),
    ("CATATATATC", "CATATC", [Variant(1, 9, "ATAT")], "2_9AT[2]"),
    ("CATATATATC", "CATATATC", [Variant(1, 9, "ATATAT")], "2_9AT[3]"),
    ("CATATATC", "CATATATATC", [Variant(1, 7, "ATATATAT")], "2_7AT[4]"),
    ("CATC", "CATATC", [Variant(1, 3, "ATAT")], "2_3dup"),
    ("CATATC", "CC", [Variant(1, 5, "")], "2_5del"),
    ("CC", "CATATC", [Variant(1, 1, "ATAT")], "1_2insAT[2]"),
    ("TCAT", "TCATCAT", [Variant(0, 4, "TCATCAT")], "2_4dup"),
    ("TCATCAT", "TCAT", [Variant(0, 7, "TCAT")], "1_6TCA[1]"),
    ("CGACTGACGTTACCGAAGTTTTTTGTACAGTCGACTGACGTTCGTCCATGATACAGAGTATGCGCAATTCC",
     "CGACTGACATTACCGAAGTTTTTTTGTACAGGGTTCTGACGATCGTCCATGGCACGGGTATGCGCGCAATTGC",
     [Variant(8, 9, "A"), Variant(18, 24, "TTTTTTT"), Variant(29, 35, "GGGTTC"), Variant(40, 41, "A"), Variant(50, 57, "GCACGG"), Variant(61, 65, "GCGCGC"), Variant(69, 70, "G")],
     "[9G>A;19_24T[7];31_34delinsGGTT;41T>A;51_57delinsGCACGG;62_65GC[3];70C>G]"),
    ("TCATCA", "CATCATCAT", [Variant(0, 0, "CA"), Variant(6, 6, "T")], "[0_1insCA;6_7insT]"),
    ("TCATCATC", "CATCATCAT", [Variant(0, 8, "CATCATCAT")], "1_8delinsCAT[3]"),
    ("CATCA", "CATCATCATCATCATCAT", [Variant(5, 5, "TCATCATCATCAT")], "5_6ins[TCA[4];T]"),
    ("CATCAT", "CATCATCATCATC", [Variant(6, 6,  "CATCATC")], "6_7ins[CAT[2];C]"),
    ("CATCAT", "CATCATCATCATCA", [Variant(6, 6,  "CATCATCA")], "6_7ins[CAT[2];CA]"),
    ("CCTTGTCTCTTATGGGATTACATCCATGG", "CCTTGGATGTAATCCCATAAGAGACATGG", [Variant(4, 25, "GGATGTAATCCCATAAGAGAC")], "6_24inv"),
    ("TACATATATATATGCATATATATATGCATATATATATATATATGCATATATATA", "TACATATATATATGCATATATATATATGCATATATATATATATGCATATATATA", [Variant(15, 43, "ATATATATATATGCATATATATATATAT")], "26_29inv"),
    ("", "", [], "="),
    ("A", "", [Variant(0, 1, "")], "1del"),
    ("", "A", [Variant(0, 0, "A")], "0_1insA"),
    ("CATCAT", "", [Variant(0, 6, "")], "1_6del"),
    ("TAA", "T", [Variant(1, 3, "")], "2_3del"),
    ("TAACCAG", "GATCC", [Variant(0, 3, "GAT"), Variant(5, 7, "")], "[1_3delinsGAT;6_7del]"),
    ("CGTGA", "CCGGATATT", [Variant(0, 1, "CC"), Variant(2, 3, ""), Variant(5, 5, "TATT")], "[1dup;3del;5_6insTATT]"),
    ("GAAG", "GGAAGCACG", [Variant(0, 1, "GG"), Variant(3, 4, "GCACG")], "[1dup;4_5insCACG]"),
    ("GAAGC", "GGAAGCACGC", [Variant(0, 1, "GG"), Variant(3, 5, "GCACGC")], "[1dup;5_6insACGC]"),
    ("CTGAAT", "TCAAATTG", [Variant(0, 3, "TCA"), Variant(6, 6, "TG")], "[1_3delinsTCA;6_7insTG]"),
    ("AGATAGCCTAACGT", "AGCCT", [Variant(0, 6, "AG"), Variant(8, 14, "T")], "[3_6del;10_14del]"),
    ("GAGTTA", "AGGTATG", [Variant(0, 1, ""), Variant(3, 4, "G"), Variant(6, 6, "TG")], "[1del;4T>G;6_7insTG]"),
    ("AGCCCCA", "AGAGGCGCC", [Variant(1, 7, "GAGGCGCC")], "3_7delinsAGGCGCC"),
    ("TCTGGAAACACTGGT", "GCGAACTAGGT", [Variant(0, 4, "GC"), Variant(6, 10, "A"), Variant(12, 12, "A")], "[1_4delinsGC;8_10del;12_13insA]"),
])
def test_extract_sequence(reference, observed, variants, hgvs):
    canonical, _ = extract_sequence(reference, observed)
    assert canonical == variants
    assert to_hgvs(canonical, reference) == hgvs


@pytest.mark.parametrize("reference, supremal, variants, hgvs", [
    ("GTGTGTTTTTTTAACAGGGA", Variant(5, 12, "TATAT"), [Variant(6, 11, "ATA")], "7_11delinsATA"),
])
def test_extract_supremal(reference, supremal, variants, hgvs):
    canonical, _ = extract_supremal(reference, supremal)
    assert canonical == variants
    assert to_hgvs(canonical, reference) == hgvs


@pytest.mark.parametrize("reference, variants, hgvs", [
    ("GTGTGTTTTTTTAACAGGGA", [Variant(6, 11, "ATA")], "7_11delinsATA"),
    ("TGCATTAGGGCAAGGGTCTTCGACTTTCCACGAAAATCGCGTCGGTTTGAC", [Variant(24, 25, "")], "25_27T[2]"),
])
def test_extract(reference, variants, hgvs):
    extracted, _ = extract(reference, variants)
    assert to_hgvs(extracted, reference) == hgvs


@pytest.mark.parametrize("reference, variants, exception, message", [
    ("", [Variant(0, 0, "")], ValueError, "empty variant"),
])
def test_to_hgvs_fail(reference, variants, exception, message):
    with pytest.raises(exception) as exc:
        print(to_hgvs(variants, reference))
    assert str(exc.value) == message
