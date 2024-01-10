"""Functions to compare variants as sequences."""


from itertools import product
from .relation import Relation
from ..lcs import LCSgraph, edit_distance


def are_equivalent(_reference, lhs, rhs):
    """Check if two variants are equivalent."""
    return lhs == rhs


def contains(reference, lhs, rhs):
    """Check if `lhs` contains `rhs`."""
    if lhs == rhs:
        return False

    lhs_distance = edit_distance(reference, lhs)
    rhs_distance = edit_distance(reference, rhs)
    distance = edit_distance(lhs, rhs)

    if lhs_distance + rhs_distance == distance:
        return False

    return lhs_distance - rhs_distance == distance


def is_contained(reference, lhs, rhs):
    """Check if `lhs` is contained in `rhs`."""
    return contains(reference, rhs, lhs)


def are_disjoint(reference, lhs, rhs):
    """Check if two variants are disjoint."""
    if lhs == rhs:
        return False

    lhs_distance = edit_distance(reference, lhs)
    rhs_distance = edit_distance(reference, rhs)
    distance = edit_distance(lhs, rhs)

    if lhs_distance + rhs_distance == distance:
        return True

    if distance == abs(lhs_distance - rhs_distance):
        return False

    lhs_graph = LCSgraph.from_sequence(reference, lhs)
    rhs_graph = LCSgraph.from_sequence(reference, rhs)

    for lhs_variant, rhs_variant in product(lhs_graph.edges(), rhs_graph.edges()):
        if not lhs_variant.is_disjoint(rhs_variant):
            return False

    return True


def have_overlap(reference, lhs, rhs):
    """Check if two variants overlap."""
    if lhs == rhs:
        return False

    lhs_distance = edit_distance(reference, lhs)
    rhs_distance = edit_distance(reference, rhs)
    distance = edit_distance(lhs, rhs)

    if distance in (lhs_distance + rhs_distance, abs(lhs_distance - rhs_distance)):
        return False

    lhs_graph = LCSgraph.from_sequence(reference, lhs)
    rhs_graph = LCSgraph.from_sequence(reference, rhs)

    for lhs_variant, rhs_variant in product(lhs_graph.edges(), rhs_graph.edges()):
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

    lhs_distance = edit_distance(reference, lhs)
    rhs_distance = edit_distance(reference, rhs)
    distance = edit_distance(lhs, rhs)

    if lhs_distance + rhs_distance == distance:
        return Relation.DISJOINT

    if lhs_distance - rhs_distance == distance:
        return Relation.CONTAINS

    if rhs_distance - lhs_distance == distance:
        return Relation.IS_CONTAINED

    lhs_graph = LCSgraph.from_sequence(reference, lhs)
    rhs_graph = LCSgraph.from_sequence(reference, rhs)

    for lhs_variant, rhs_variant in product(lhs_graph.edges(), rhs_graph.edges()):
        if not lhs_variant.is_disjoint(rhs_variant):
            return Relation.OVERLAP

    return Relation.DISJOINT
