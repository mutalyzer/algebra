"""A Boolean Algebra for Genetic Variants."""


from .relations import Relation
from .relations.variant_based import compare
from .variants import Variant


__all__ = [
    "Relation",
    "Variant",
    "compare",
]
