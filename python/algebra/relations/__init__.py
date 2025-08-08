"""Calculate the relation between two variants using the definitions
from [1]_. Variants are given as alleles, i.e., a list of variants
applied to some reference sequence. Both variants are subject to the same
reference sequence.

References
----------
[1] J.K. Vis, M.A. Santcroos, W.A. Kosters and J.F.J. Laros.
"A Boolean Algebra for Genetic Variants." In: Bioinformatics (2023).
"""


from algebra_ext import Relation, compare


__all__ = [
    "Relation",
    "compare",
]
