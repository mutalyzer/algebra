"""Extract the canonical variant."""


from .extractor import extract, extract_supremal, to_hgvs
from .repeats import extract_repeats


__all__ = [
    "extract",
    "extract_repeats",
    "extract_supremal",
    "to_hgvs",
]
