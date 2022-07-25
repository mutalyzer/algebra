"""Relation class and functions to compare variants.

Calculate the relation between two variants using the definitions
from [1]_. Variants are given as observed sequences (alleles), i.e., a
list of variants applied to some reference sequence. Both variants
subject to the same reference sequence.

References
----------
[1] J.K. Vis, M.A. Santcroos, W.A. Kosters and J.F.J. Laros.
"A Boolean Algebra for Genetic Variants".
In: arXiv preprint 2112.14494 (2021).
"""


from enum import Enum


class Relation(Enum):
    """Relation enum."""
    EQUIVALENT = "equivalent"
    CONTAINS = "contains"
    IS_CONTAINED = "is_contained"
    OVERLAP = "overlap"
    DISJOINT = "disjoint"
