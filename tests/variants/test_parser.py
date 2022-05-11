import pytest
from algebra.variants import Parser, Variant


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
    ("4dupT", [Variant(4, 4, "T")]),
    ("0_1insT", [Variant(0, 0, "T")]),
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
    ("3=", ValueError, "expected '>' at 1"),
    ("[3del", ValueError, "expected ']' at 5"),
    ("[3del;", ValueError, "unexpected end of expression"),
    ("3del;", ValueError, "expected end of expression at 4"),
    ("4dupTT", ValueError, "inconsistent duplicated length at 6"),
    ("4dup", NotImplementedError, ""),
    ("4inv", NotImplementedError, ""),
])
def test_hgvs_parser_fail(expression, exception, message):
    with pytest.raises(exception) as exc:
        Parser(expression).hgvs()
    assert str(exc.value) == message


@pytest.mark.parametrize("expression, variants", [
    ("AAA:0:0:", [Variant(0, 0)]),
    ("AAA:0:3:TTT", [Variant(0, 3, "TTT")]),
    ("AAA:0:AAA:TTT", [Variant(0, 3, "TTT")]),
])
def test_spdi_parser(expression, variants):
    assert Parser(expression).spdi() == variants


@pytest.mark.parametrize("expression, exception, message", [
    ("", ValueError, "unexpected end of expression"),
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
