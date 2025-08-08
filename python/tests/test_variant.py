import pytest
from algebra.variants import Variant


@pytest.mark.parametrize("args, expected", [
    ((0, 0, ""), {"start": 0, "end": 0, "sequence": ""}),
    ((0, 0, "A"), {"start": 0, "end": 0, "sequence": "A"}),
    ((0, 1, "A"), {"start": 0, "end": 1, "sequence": "A"}),
])
def test_variant(args, expected):
    variant = Variant(*args)
    assert all((variant.start == expected["start"],
        variant.end == expected["end"],
        variant.sequence == expected["sequence"]))


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


@pytest.mark.parametrize("args, expected", [
    ((0, 0, ""), "0:0/"),
    ((0, 0, "A"), "0:0/A"),
    ((0, 1, "A"), "0:1/A"),
])
def test_variant_repr(args, expected):
    assert str(Variant(*args)) == expected


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
