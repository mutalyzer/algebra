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
from itertools import product
from .lcs import edit, edit_distance_only, lcs_graph


class Relation(Enum):
    """Relation enum."""
    EQUIVALENT = "equivalent"
    CONTAINS = "contains"
    IS_CONTAINED = "is_contained"
    OVERLAP = "overlap"
    DISJOINT = "disjoint"


def are_equivalent(_reference, lhs, rhs):
    """Check if two variants are equivalent."""
    return lhs == rhs


def contains(reference, lhs, rhs):
    """Check if `lhs` contains `rhs`."""
    if lhs == rhs:
        return False

    lhs_distance = edit_distance_only(reference, lhs)
    rhs_distance = edit_distance_only(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if lhs_distance + rhs_distance == distance:
        return False

    return lhs_distance - rhs_distance == distance


def is_contained(reference, lhs, rhs):
    """Check if `lhs` is contained in `rhs`."""
    if lhs == rhs:
        return False

    lhs_distance = edit_distance_only(reference, lhs)
    rhs_distance = edit_distance_only(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if lhs_distance + rhs_distance == distance:
        return False

    return rhs_distance - lhs_distance == distance


def are_disjoint(reference, lhs, rhs):
    """Check if two variants are disjoint."""
    if lhs == rhs:
        return False

    lhs_distance, lhs_lcs_nodes = edit(reference, lhs)
    rhs_distance, rhs_lcs_nodes = edit(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if lhs_distance + rhs_distance == distance:
        return True

    if distance == abs(lhs_distance - rhs_distance):
        return False

    _, lhs_edges = lcs_graph(reference, lhs, lhs_lcs_nodes)
    _, rhs_edges = lcs_graph(reference, rhs, rhs_lcs_nodes)

    for lhs_variant, rhs_variant in product(lhs_edges, rhs_edges):
        if not lhs_variant.is_disjoint(rhs_variant):
            return False

    return True


def have_overlap(reference, lhs, rhs):
    """Check if two variants overlap."""
    if lhs == rhs:
        return False

    lhs_distance, lhs_lcs_nodes = edit(reference, lhs)
    rhs_distance, rhs_lcs_nodes = edit(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if distance in (lhs_distance + rhs_distance, abs(lhs_distance - rhs_distance)):
        return False

    _, lhs_edges = lcs_graph(reference, lhs, lhs_lcs_nodes)
    _, rhs_edges = lcs_graph(reference, rhs, rhs_lcs_nodes)

    for lhs_variant, rhs_variant in product(lhs_edges, rhs_edges):
        if not lhs_variant.is_disjoint(rhs_variant):
            return True

    return False


def compare(reference, lhs, rhs):
    """Compare two variants.

    Parameters
    ----------
    reference : str
        The reference sequence.
    lhs : str
        The variant (as observed sequence) on the left-hand side.
    rhs : str
        The variant on the right-hand side.

    Returns
    -------
    `Relation`
        The relation between the two variants.
    """

    if lhs == rhs:
        return Relation.EQUIVALENT

    lhs_distance, lhs_lcs_nodes = edit(reference, lhs)
    rhs_distance, rhs_lcs_nodes = edit(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if lhs_distance + rhs_distance == distance:
        return Relation.DISJOINT

    if lhs_distance - rhs_distance == distance:
        return Relation.CONTAINS

    if rhs_distance - lhs_distance == distance:
        return Relation.IS_CONTAINED

    _, lhs_edges = lcs_graph(reference, lhs, lhs_lcs_nodes)
    _, rhs_edges = lcs_graph(reference, rhs, rhs_lcs_nodes)

    for lhs_variant, rhs_variant in product(lhs_edges, rhs_edges):
        if not lhs_variant.is_disjoint(rhs_variant):
            return Relation.OVERLAP

    return Relation.DISJOINT
