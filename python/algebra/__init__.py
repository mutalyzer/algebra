"""A Boolean Algebra for Genetic Variants."""


from .lcs import LCSgraph
from .relations import Relations, compare
from .variants import Variant


__all__ = [
    "LCSgraph",
    "Relations",
    "Variant",
    "compare",
]
