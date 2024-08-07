import pytest
from algebra import Variant
from algebra.variants import patch, to_hgvs


@pytest.mark.parametrize("args, expected", [
    ((0, 0, ""), Variant(0, 0, "")),
    ((10, 12, "T"), Variant(10, 12, "T")),
])
def test_variant(args, expected):
    assert Variant(*args) == expected


@pytest.mark.parametrize("args, exception, message", [
    (("1", "2", 3), TypeError, "start must be an integer"),
    ((1, "2", 3), TypeError, "end must be an integer"),
    ((1, 2, 3), TypeError, "sequence must be a string"),
    ((2, 1, "A"), ValueError, "start must not be after end"),
    ((-1, 0, ""), ValueError, "start must be greater or equal to 0"),
])
def test_variant_fail(args, exception, message):
    with pytest.raises(exception) as exc:
        Variant.create_safe(*args)
    assert str(exc.value) == message


@pytest.mark.parametrize("variant", [
    (Variant(0, 1, "")),
    (Variant(0, 0, "A")),
])
def test_variant_bool_true(variant):
    assert variant


@pytest.mark.parametrize("variant", [
    (Variant(0, 0, "")),
])
def test_variant_bool_false(variant):
    assert not variant


@pytest.mark.parametrize("lhs, rhs", [
    (Variant(0, 0, "A"), Variant(0, 0, "A")),
])
def test_variant_equal_true(lhs, rhs):
    assert lhs == rhs


@pytest.mark.parametrize("lhs, rhs", [
    (Variant(0, 0, "T"), Variant(0, 0, "A")),
    (Variant(0, 1, ""), Variant(0, 0, "")),
])
def test_variant_equal_false(lhs, rhs):
    assert not lhs == rhs


@pytest.mark.parametrize("duplicates, unique", [
    ([Variant(0, 0, ""), Variant(0, 0, "")], {Variant(0, 0, "")}),
])
def test_variant_hash(duplicates, unique):
    assert set(duplicates) == unique


@pytest.mark.parametrize("variant, distance", [
    (Variant(0, 1, "A"), 2),
])
def test_variant_len(variant, distance):
    assert len(variant) == distance


@pytest.mark.parametrize("variants, expected", [
    ([Variant(0, 0, "C"), Variant(0, 1, "C")], [Variant(0, 0, "C"), Variant(0, 1, "C")]),
    ([Variant(0, 1, "C"), Variant(0, 0, "C")], [Variant(0, 0, "C"), Variant(0, 1, "C")]),
    ([Variant(3, 4, ""), Variant(1, 2, "")], [Variant(1, 2, ""), Variant(3, 4, "")]),
    ([Variant(2, 4, "T"), Variant(2, 2, "G")], [Variant(2, 2, "G"), Variant(2, 4, "T")]),
    ([Variant(4, 4, "G"), Variant(2, 4, "T")], [Variant(2, 4, "T"), Variant(4, 4, "G")]),
])
def test_variant_sort(variants, expected):
    assert sorted(variants) == expected


@pytest.mark.parametrize("variants", [
    ([Variant(1, 3, "C"), Variant(0, 2, "")]),
    ([Variant(4, 4, "C"), Variant(4, 4, "C")]),
    ([Variant(2, 4, "T"), Variant(3, 3, "G")]),
    ([Variant(0, 0, ""), Variant(0, 0, "")]),
])
def test_variant_sort_fail(variants):
    with pytest.raises(ValueError) as exc:
        sorted(variants)
    assert str(exc.value) == "unorderable variants"


@pytest.mark.parametrize("variant, string", [
    (Variant(0, 0, ""), "0:0/"),
    (Variant(0, 0, "TTT"), "0:0/TTT"),
    (Variant(0, 1, "T"), "0:1/T"),
])
def test_variant_string(variant, string):
    assert str(variant) == string


@pytest.mark.parametrize("lhs, rhs, expected", [
    (Variant(0, 1, ""), Variant(4, 5, ""), True),
    (Variant(0, 0, "T"), Variant(4, 5, "T"), True),
    (Variant(0, 0, "T"), Variant(0, 1, "T"), False),
    (Variant(0, 5, ""), Variant(2, 3, ""), False),
    (Variant(0, 0, "C"), Variant(0, 1, "T"), True),
    (Variant(2, 3, "T"), Variant(2, 2, "T"), False),
    (Variant(2, 3, "T"), Variant(3, 3, "T"), False),
    (Variant(2, 3, "T"), Variant(2, 2, "C"), True),
    (Variant(2, 3, "T"), Variant(3, 3, "C"), True),
    (Variant(2, 2, "T"), Variant(2, 2, "T"), False),
    (Variant(2, 3, ""), Variant(1, 2, ""), True),
    (Variant(2, 3, ""), Variant(3, 4, ""), True),
    (Variant(2, 3, "T"), Variant(2, 2, "C"), True),
    (Variant(2, 3, "T"), Variant(3, 3, "C"), True),
    (Variant(2, 3, "CT"), Variant(3, 3, "C"), False),
])
def test_variant_is_disjoint(lhs, rhs, expected):
    assert lhs.is_disjoint(rhs) == rhs.is_disjoint(lhs) == expected


@pytest.mark.parametrize("variant, hgvs", [
    (Variant(0, 0, ""), "="),
    (Variant(5, 5, ""), "="),
    (Variant(2, 3, ""), "3del"),
    (Variant(2, 4, ""), "3_4del"),
    (Variant(2, 3, "AA"), "3delinsAA"),
    (Variant(2, 4, "TT"), "3_4delinsTT"),
    (Variant(3, 3, "TTT"), "3_4insTTT"),
    (Variant(4, 5, "G"), "5>G"),
])
def test_variant_to_hgvs(variant, hgvs):
    assert variant.to_hgvs() == hgvs


@pytest.mark.parametrize("reference, variant, hgvs", [
    ("AAAAAA", Variant(4, 5, "G"), "5A>G"),
])
def test_variant_to_hgvs_reference(reference, variant, hgvs):
    assert variant.to_hgvs(reference) == hgvs


@pytest.mark.parametrize("reference, variant, hgvs", [
    ("AAAAAA", Variant(3, 5, ""), "4_5delAA"),
    ("AAAAAA", Variant(3, 5, "T"), "4_5delAAinsT"),
])
def test_variant_to_hgvs_deletions(reference, variant, hgvs):
    assert variant.to_hgvs(reference, only_substitutions=False) == hgvs


@pytest.mark.parametrize("reference, variant, spdi", [
    ("AAA", Variant(10, 12, "TT"), "AAA:10:2:TT"),
    ("AAA", Variant(10, 10, ""), "AAA:10:0:"),
])
def test_variant_to_spdi(reference, variant, spdi):
    assert variant.to_spdi(reference) == spdi


@pytest.mark.parametrize("variant, atomics", [
    (Variant(3, 6, "XYZ"), [
        [Variant(3, 3, "XYZ"), Variant(3, 4, ""), Variant(4, 5, ""), Variant(5, 6, "")],
        [Variant(3, 3, "XY"), Variant(3, 4, ""), Variant(4, 4, "Z"), Variant(4, 5, ""), Variant(5, 6, "")],
        [Variant(3, 3, "XY"), Variant(3, 4, ""), Variant(4, 5, ""), Variant(5, 5, "Z"), Variant(5, 6, "")],
        [Variant(3, 3, "XY"), Variant(3, 4, ""), Variant(4, 5, ""), Variant(5, 6, ""), Variant(6, 6, "Z")],
        [Variant(3, 3, "X"), Variant(3, 4, ""), Variant(4, 4, "YZ"), Variant(4, 5, ""), Variant(5, 6, "")],
        [Variant(3, 3, "X"), Variant(3, 4, ""), Variant(4, 4, "Y"), Variant(4, 5, ""), Variant(5, 5, "Z"), Variant(5, 6, "")],
        [Variant(3, 3, "X"), Variant(3, 4, ""), Variant(4, 4, "Y"), Variant(4, 5, ""), Variant(5, 6, ""), Variant(6, 6, "Z")],
        [Variant(3, 3, "X"), Variant(3, 4, ""), Variant(4, 5, ""), Variant(5, 5, "YZ"), Variant(5, 6, "")],
        [Variant(3, 3, "X"), Variant(3, 4, ""), Variant(4, 5, ""), Variant(5, 5, "Y"), Variant(5, 6, ""), Variant(6, 6, "Z")],
        [Variant(3, 3, "X"), Variant(3, 4, ""), Variant(4, 5, ""), Variant(5, 6, ""), Variant(6, 6, "YZ")],
        [Variant(3, 4, ""), Variant(4, 4, "XYZ"), Variant(4, 5, ""), Variant(5, 6, "")],
        [Variant(3, 4, ""), Variant(4, 4, "XY"), Variant(4, 5, ""), Variant(5, 5, "Z"), Variant(5, 6, "")],
        [Variant(3, 4, ""), Variant(4, 4, "XY"), Variant(4, 5, ""), Variant(5, 6, ""), Variant(6, 6, "Z")],
        [Variant(3, 4, ""), Variant(4, 4, "X"), Variant(4, 5, ""), Variant(5, 5, "YZ"), Variant(5, 6, "")],
        [Variant(3, 4, ""), Variant(4, 4, "X"), Variant(4, 5, ""), Variant(5, 5, "Y"), Variant(5, 6, ""), Variant(6, 6, "Z")],
        [Variant(3, 4, ""), Variant(4, 4, "X"), Variant(4, 5, ""), Variant(5, 6, ""), Variant(6, 6, "YZ")],
        [Variant(3, 4, ""), Variant(4, 5, ""), Variant(5, 5, "XYZ"), Variant(5, 6, "")],
        [Variant(3, 4, ""), Variant(4, 5, ""), Variant(5, 5, "XY"), Variant(5, 6, ""), Variant(6, 6, "Z")],
        [Variant(3, 4, ""), Variant(4, 5, ""), Variant(5, 5, "X"), Variant(5, 6, ""), Variant(6, 6, "YZ")],
        [Variant(3, 4, ""), Variant(4, 5, ""), Variant(5, 6, ""), Variant(6, 6, "XYZ")]
    ]),
])
def test_variant_atomics(variant, atomics):
    assert list(variant.atomics()) == atomics


@pytest.mark.parametrize("variant, end, expected", [
    (Variant(7, 8, "G"), 10, Variant(1, 2, "C")),
])
def test_variant_reverse_complement(variant, end, expected):
    assert variant.reverse_complement(end) == expected


@pytest.mark.parametrize("reference, variants, observed", [
    ("ACCTGC", [Variant(1, 4, "CCC")], "ACCCGC"),
    ("CAAAAC", [Variant(1, 5, "AAAAAAAA")], "CAAAAAAAAC"),
])
def test_patch(reference, variants, observed):
    assert patch(reference, variants) == observed


@pytest.mark.parametrize("reference, variants, hgvs", [
    ("AAA", [], "AAA:g.="),
    ("ACCTGC", [Variant(1, 4, "CCC")], "ACCTGC:g.2_4delinsCCC"),
    ("ACCTGC", [Variant(3, 4, "C"), Variant(4, 5, "T")], "ACCTGC:g.[4T>C;5G>T]"),
])
def test_to_hgvs(reference, variants, hgvs):
    assert to_hgvs(variants, reference, sequence_prefix=True) == hgvs
