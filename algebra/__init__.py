"""A Boolean Algebra for Genetic Variants."""


from .relations import Relation
from .relations.variant_based import are_disjoint, are_equivalent, compare, contains, have_overlap, is_contained
from .variants import Variant


__all__ = [
    "Relation",
    "Variant",
    "are_disjoint",
    "are_equivalent",
    "compare",
    "contains",
    "have_overlap",
    "is_contained",
]
