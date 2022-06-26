import pytest
from algebra.variants import Parser, Variant, reverse_complement


@pytest.mark.parametrize("expression, variants", [
    ("=", []),
    ("3del", [Variant(2, 3)]),
    ("3delA", [Variant(2, 3)]),
    ("3_4del", [Variant(2, 4)]),
    ("3_4delTT", [Variant(2, 4)]),
    ("3delinsA", [Variant(2, 3, "A")]),
    ("3delAinsA", [Variant(2, 3, "A")]),
    ("3_4delinsTT", [Variant(2, 4, "TT")]),
    ("3_4delTTinsTT", [Variant(2, 4, "TT")]),
    ("3_4insTTT", [Variant(3, 3, "TTT")]),
    ("3A>T", [Variant(2, 3, "T")]),
    ("5>G", [Variant(4, 5, "G")]),
    ("[3del]", [Variant(2, 3)]),
    ("[3del;3_4insT]", [Variant(2, 3), Variant(3, 3, "T")]),
    ("3_4A[4]", [Variant(2, 4, "AAAA")]),
    ("0_1insT", [Variant(0, 0, "T")]),
    ("3=", []),
    ("[1=;2=;3=]", []),
])
def test_hgvs_parser(expression, variants):
    assert Parser(expression).hgvs() == variants


@pytest.mark.parametrize("expression, exception, message", [
    (None, TypeError, "expression must be a string"),
    ("", ValueError, "unexpected end of expression"),
    ("0del", ValueError, "start must be greater or equal to 0"),
    ("del", ValueError, "expected digit at 0"),
    ("3_3del", ValueError, "invalid range at 3"),
    ("4_3del", ValueError, "invalid range at 3"),
    ("3delAA", ValueError, "inconsistent deleted length at 6"),
    ("3_4delA", ValueError, "inconsistent deleted length at 7"),
    ("3insA", ValueError, "invalid inserted range at 4"),
    ("3_5insA", ValueError, "invalid inserted range at 6"),
    ("3_4ins", ValueError, "unexpected end of expression"),
    ("10_12", ValueError, "unexpected end of expression"),
    ("10_12[", ValueError, "expected nucleotide at 5"),
    ("10_12A", ValueError, "expected '[' at 6"),
    ("10_12A[", ValueError, "unexpected end of expression"),
    ("10_12A[A", ValueError, "expected digit at 7"),
    ("10_12A[1", ValueError, "expected ']' at 8"),
    ("123", ValueError, "expected '>' at 3"),
    ("=3", ValueError, "expected end of expression at 1"),
    ("[3del", ValueError, "expected ']' at 5"),
    ("[3del;", ValueError, "unexpected end of expression"),
    ("3del;", ValueError, "expected end of expression at 4"),
    ("4dup", NotImplementedError, "duplications without the reference context are not supported"),
    ("4inv", NotImplementedError, "inversions without the reference context are not supported"),
    ("6875TTTCGCCCC[3", NotImplementedError, "dbSNP repeats without reference context are not supported"),
])
def test_hgvs_parser_fail(expression, exception, message):
    with pytest.raises(exception) as exc:
        Parser(expression).hgvs()
    assert str(exc.value) == message


@pytest.mark.parametrize("expression, variants", [
    ("NG_008376.4:g.=", []),
    ("NG_008376.4:g.3del", [Variant(2, 3)]),
])
def test_hgvs_parser_with_prefix(expression, variants):
    assert Parser(expression).hgvs(skip_reference=True) == variants


@pytest.mark.parametrize("expression, exception, message", [
    ("3del", ValueError, "expected 'g.' at 4"),
    ("g.3del", ValueError, "expected 'g.' at 6"),
    ("NG_008376.4:", ValueError, "expected 'g.' at 12"),
    ("NG_008376.4:3del", ValueError, "expected 'g.' at 12"),
    ("NG_008376.4:g.", ValueError, "unexpected end of expression"),
])
def test_hgvs_parser_with_prefix_fail(expression, exception, message):
    with pytest.raises(exception) as exc:
        Parser(expression).hgvs(skip_reference=True)
    assert str(exc.value) == message


@pytest.mark.parametrize("reference, expression, variants", [
    ("ACCGGGTTTT", "1inv", [Variant(0, 1, "T")]),
    ("ACCGGGTTTT", "1_10inv", [Variant(0, 10, "AAAACCCGGT")]),
    ("ACCGGGTTTT", "1dup", [Variant(1, 1, "A")]),
    ("ACCGGGTTTT", "1_2dup", [Variant(2, 2, "AC")]),
    ("TTGAGAGAGATT", "3GA[3]", [Variant(2, 10, "GAGAGA")]),
    ("AAA", "1delA", [Variant(0, 1)]),
])
def test_hgvs_parser_with_reference(reference, expression, variants):
    assert Parser(expression).hgvs(reference) == variants


@pytest.mark.parametrize("reference, expression, exception, message", [
    ("ACCGGGTTTT", "1_11inv", ValueError, "invalid range in reference"),
    ("ACCGGGTTTT", "11dup", ValueError, "invalid range in reference"),
    ("ACCGGGTTTT", "0_1dup", ValueError, "invalid range in reference"),
    ("TTGAGAGAGATT", "3GA[3", ValueError, "expected ']' at 5"),
    ("TTGAGAGAGATT", "3AG[3]", ValueError, "'AG' not found in reference at 2"),
    ("AAAAAA", "4dupTT", ValueError, "inconsistent duplicated length at 6"),
    ("AAAAAA", "4dupT", ValueError, "'T' not found in reference at 3"),
    ("AAA", "1delT", ValueError, "'T' not found in reference at 0"),
])
def test_hgvs_parser_with_reference_fail(reference, expression, exception, message):
    with pytest.raises(exception) as exc:
        Parser(expression).hgvs(reference)
    assert str(exc.value) == message


@pytest.mark.parametrize("expression, variants", [
    ("AAA:0:0:", [Variant(0, 0)]),
    ("AAA:0:3:TTT", [Variant(0, 3, "TTT")]),
    ("AAA:0:AAA:TTT", [Variant(0, 3, "TTT")]),
    (":1:TT:G", [Variant(1, 3, "G")]),
])
def test_spdi_parser(expression, variants):
    assert Parser(expression).spdi() == variants


@pytest.mark.parametrize("expression, exception, message", [
    ("", ValueError, "expected ':' at 0"),
    ("A", ValueError, "expected ':' at 1"),
    ("A:", ValueError, "unexpected end of expression"),
    ("A:0", ValueError, "expected ':' at 3"),
    ("A:0:", ValueError, "expected ':' at 4"),
    ("A:0:0", ValueError, "expected ':' at 5"),
    ("A:0:0::", ValueError, "expected end of expression at 6"),
])
def test_spdi_parser_fail(expression, exception, message):
    with pytest.raises(exception) as exc:
        Parser(expression).spdi()
    assert str(exc.value) == message


@pytest.mark.parametrize("sequence, expected", [
    ("", ""),
    ("A", "T"),
    ("ACGT", "ACGT"),
    ("ACCGGGTTTT", "AAAACCCGGT"),
])
def test_reverse_complement(sequence, expected):
    assert reverse_complement(sequence) == expected
