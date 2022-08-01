"""Relation class (enum)."""


from enum import Enum


class Relation(Enum):
    """Relation enum."""
    EQUIVALENT = "equivalent"
    CONTAINS = "contains"
    IS_CONTAINED = "is_contained"
    OVERLAP = "overlap"
    DISJOINT = "disjoint"
