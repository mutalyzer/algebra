import pytest
from algebra import Variant
from algebra.variants import parse_hgvs, parse_spdi, reverse_complement


@pytest.mark.parametrize("expression, variants", [
    ("=", []),
    ("3del", [Variant(2, 3, "")]),
    ("3delA", [Variant(2, 3, "")]),
    ("3_4del", [Variant(2, 4, "")]),
    ("3_3del", [Variant(2, 3, "")]),
    ("3_4delTT", [Variant(2, 4, "")]),
    ("3delinsA", [Variant(2, 3, "A")]),
    ("3delAinsA", [Variant(2, 3, "A")]),
    ("3_4delinsTT", [Variant(2, 4, "TT")]),
    ("3_4delTTinsTT", [Variant(2, 4, "TT")]),
    ("3_4insTTT", [Variant(3, 3, "TTT")]),
    ("3A>T", [Variant(2, 3, "T")]),
    ("5>G", [Variant(4, 5, "G")]),
    ("[3del]", [Variant(2, 3, "")]),
    ("[3del;3_4insT]", [Variant(2, 3, ""), Variant(3, 3, "T")]),
    ("3_4A[4]", [Variant(2, 4, "AAAA")]),
    ("0_1insT", [Variant(0, 0, "T")]),
    ("3=", []),
    ("[1=;2=;3=]", []),
    ("g.3del", [Variant(2, 3, "")]),
    ("NG_008376.4:g.=", []),
    ("NG_008376.4:g.3del", [Variant(2, 3, "")]),
    ("NG_008376.4:3del", [Variant(2, 3, "")]),
    ("3_4invAA", [Variant(2, 4, "AA")]),
    ("6_7ins[CAT[2];C]", [Variant(6, 6, "CATCATC")]),
    ("6_7ins[A]", [Variant(6, 6, "A")]),
    ("6delins[A;A]", [Variant(5, 6, "AA")]),
    ("6delins[A[0]]", [Variant(5, 6, "")]),
])
def test_hgvs_parser(expression, variants):
    assert parse_hgvs(expression) == variants


@pytest.mark.parametrize("expression, exception, message", [
    (None, TypeError, "expression must be a string"),
    ("", ValueError, "unexpected end of expression"),
    ("0del", ValueError, "start must be greater or equal to 0"),
    ("del", ValueError, "expected digit at 1"),
    ("4_3del", ValueError, "invalid range at 3"),
    ("3delAA", ValueError, "inconsistent deleted length at 6"),
    ("3_4delA", ValueError, "inconsistent deleted length at 7"),
    ("3insA", ValueError, "invalid inserted range at 4"),
    ("3_3insT", ValueError, "invalid inserted range at 6"),
    ("3_5insA", ValueError, "invalid inserted range at 6"),
    ("3_4ins", ValueError, "unexpected end of expression"),
    ("10_12", NotImplementedError, "unsupported variant at 6"),
    ("10_12[", ValueError, "unexpected end of expression"),
    ("10_12A", NotImplementedError, "unsupported variant at 6"),
    ("10_12A[", ValueError, "unexpected end of expression"),
    ("10_12A[A", ValueError, "expected digit at 8"),
    ("10_12A[1", ValueError, "unexpected end of expression"),
    ("123", NotImplementedError, "unsupported variant at 4"),
    ("=3", ValueError, "expected end of expression at 2"),
    ("[3del", ValueError, "unexpected end of expression"),
    ("[3del;", ValueError, "unexpected end of expression"),
    ("3del;", ValueError, "expected end of expression at 5"),
    ("4dup", NotImplementedError, "duplication without reference context at 2"),
    ("4inv", NotImplementedError, "inversion without reference context at 2"),
    ("6875TTTCGCCCC[3]", NotImplementedError, "NCBI style repeat without reference context at 5"),
    ("NG_008376.4:", ValueError, "unexpected end of expression"),
    ("NG_008376.4:g.", ValueError, "unexpected end of expression"),
    ("3_4invA", ValueError, "inconsistent inversion length at 4"),
    ("3_4A>T", ValueError, "inconstistent deletion length at 4"),
    ("[3del];", ValueError, "expected end of expression at 7"),
    ("6_7ins[", ValueError, "unexpected end of expression"),
    ("6_7ins[A", ValueError, "unexpected end of expression"),
    ("6_7ins[;]", ValueError, "expected nucleotide at 8"),
    ("6_7ins[A;]", ValueError, "expected nucleotide at 10"),
    ("6_7ins[A;A", ValueError, "unexpected end of expression"),
    ("6_7ins[A;A[", ValueError, "unexpected end of expression"),
    ("6_7ins[A;A[0", ValueError, "unexpected end of expression"),
    ("6_7ins[A;A[0]", ValueError, "unexpected end of expression"),
])
def test_hgvs_parser_fail(expression, exception, message):
    with pytest.raises(exception) as exc:
        parse_hgvs(expression)
    assert str(exc.value) == message


@pytest.mark.parametrize("reference, expression, variants", [
    ("ACCGGGTTTT", "1inv", [Variant(0, 1, "T")]),
    ("ACCGGGTTTT", "1_10inv", [Variant(0, 10, "AAAACCCGGT")]),
    ("ACCGGGTTTT", "1dup", [Variant(0, 1, "AA")]),
    ("ACCGGGTTTT", "1_2dup", [Variant(0, 2, "ACAC")]),
    ("TTGAGAGAGATT", "3GA[3]", [Variant(2, 10, "GAGAGA")]),
    ("AAA", "1delA", [Variant(0, 1, "")]),
    ("GGGG", "2_3invCC", [Variant(1, 3, "CC")]),
    ("CAAAAC", "2_5A[8]", [Variant(1, 5, "AAAAAAAA")]),
])
def test_hgvs_parser_with_reference(reference, expression, variants):
    assert parse_hgvs(expression, reference=reference) == variants


@pytest.mark.parametrize("reference, expression, exception, message", [
    ("ACCGGGTTTT", "1_11inv", ValueError, "invalid range in reference"),
    ("ACCGGGTTTT", "11dup", ValueError, "invalid range in reference"),
    ("ACCGGGTTTT", "0_1dup", ValueError, "start must be greater or equal to 0"),
    ("TTGAGAGAGATT", "3GA[3", ValueError, "unexpected end of expression"),
    ("TTGAGAGAGATT", "3AG[3]", ValueError, "'AG' not found in reference at 2"),
    ("AAAAAA", "4dupTT", ValueError, "inconsistent duplicated length at 6"),
    ("AAAAAA", "4dupT", ValueError, "'T' not found in reference at 3"),
    ("AAA", "1delT", ValueError, "'T' not found in reference at 0"),
    ("GTTG", "2_3invCC", ValueError, "'CC' not found in reference at 1"),
    ("GGGG", "3A>T", ValueError, "'A' not found in reference at 2"),
])
def test_hgvs_parser_with_reference_fail(reference, expression, exception, message):
    with pytest.raises(exception) as exc:
        parse_hgvs(expression, reference=reference)
    assert str(exc.value) == message


@pytest.mark.parametrize("expression, variants", [
    ("AAA:0:0:", [Variant(0, 0, "")]),
    ("AAA:0:3:TTT", [Variant(0, 3, "TTT")]),
    ("AAA:0:AAA:TTT", [Variant(0, 3, "TTT")]),
    (":1:TT:G", [Variant(1, 3, "G")]),
])
def test_spdi_parser(expression, variants):
    assert parse_spdi(expression) == variants


@pytest.mark.parametrize("expression, exception, message", [
    (None, TypeError, "expression must be a string"),
    ("", ValueError, "not enough values to unpack (expected 4, got 1)"),
    ("A", ValueError, "not enough values to unpack (expected 4, got 1)"),
    ("A:", ValueError, "not enough values to unpack (expected 4, got 2)"),
    ("A:0", ValueError, "not enough values to unpack (expected 4, got 2)"),
    ("A:0:", ValueError, "not enough values to unpack (expected 4, got 3)"),
    ("A:0:0", ValueError, "not enough values to unpack (expected 4, got 3)"),
    ("A:0:0::", ValueError, "too many values to unpack (expected 4)"),
])
def test_spdi_parser_fail(expression, exception, message):
    with pytest.raises(exception) as exc:
        parse_spdi(expression)
    assert str(exc.value) == message


@pytest.mark.parametrize("sequence, expected", [
    ("", ""),
    ("A", "T"),
    ("ACGT", "ACGT"),
    ("ACCGGGTTTT", "AAAACCCGGT"),
])
def test_reverse_complement(sequence, expected):
    assert reverse_complement(sequence) == expected
