import pytest
from algebra.variants import Variant


@pytest.mark.parametrize("start, end, sequence", [
    (0, 0, ""),
    (0, 0, "A"),
    (0, 1, "A"),
])
def test_variant(start, end, sequence):
    variant = Variant(start, end, sequence)
    assert all((variant.start == start,
        variant.end == end,
        variant.sequence == sequence))


@pytest.mark.parametrize("args, expected", [
    ((), "Variant() missing required argument 'start' (pos 1)"),
    ((0,), "Variant() missing required argument 'end' (pos 2)"),
    ((0, 0), "Variant() missing required argument 'sequence' (pos 3)"),
    (("0", 0, ""), "'str' object cannot be interpreted as an integer"),
    ((0, "0", ""), "'str' object cannot be interpreted as an integer"),
    ((0, 0, 0), "a bytes-like object is required, not 'int'"),
])
def test_variant_TypeError(args, expected):
    with pytest.raises(TypeError) as exc:
        Variant(*args)
    assert str(exc.value) == expected


@pytest.mark.parametrize("variant, expected", [
    (Variant(0, 0, ""), "0:0/"),
    (Variant(0, 0, "A"), "0:0/A"),
    (Variant(0, 1, "A"), "0:1/A"),
])
def test_variant_repr(variant, expected):
    assert str(variant) == expected


def test_variant_immutable():
    variant = Variant(sequence="", start=0, end=0)

    with pytest.raises(AttributeError) as exc:
        variant.start = 1
    assert str(exc.value) == "readonly attribute"

    with pytest.raises(AttributeError) as exc:
        variant.end = 1
    assert str(exc.value) == "readonly attribute"

    with pytest.raises(AttributeError) as exc:
        variant.sequence = "A"
    assert str(exc.value) == "readonly attribute"


def test_variant_compare():
    assert Variant(0, 0, "") == Variant(0, 0, "")
    assert Variant(0, 0, "A") == Variant(0, 0, "A")
    assert Variant(0, 0, "") != Variant(0, 1, "")
    assert Variant(0, 0, "") != Variant(1, 0, "")
    assert Variant(0, 0, "") != Variant(0, 0, "A")


def test_variant_from_slice():
    reference = "GTTCGCGGGGAAAGGAAAAAAGCCGCCGGGCAGGAAA"
    variants = [Variant(1, 5, reference[7:34])]
    del reference
    assert variants == [Variant(1, 5, "GGGAAAGGAAAAAAGCCGCCGGGCAGG")]
