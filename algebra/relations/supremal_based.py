"""Functions to compare supremal variants."""


from ..lcs import LCSgraph
from ..relations import Relation
from .graph_based import (are_disjoint as graph_based_are_disjoint,
                          compare as graph_based_compare,
                          have_overlap as graph_based_have_overlap)
from .sequence_based import contains as sequence_based_contains


def are_equivalent(_reference, lhs, rhs):
    """Check if two variants are equivalent."""
    return lhs == rhs


def contains(reference, lhs, rhs):
    """Check if `lhs` contains `rhs`."""
    if lhs == rhs or not lhs or not rhs or lhs.is_disjoint(rhs):
        return False

    start = min(lhs.start, rhs.start)
    end = max(lhs.end, rhs.end)
    lhs_observed = "".join([reference[start:lhs.start],
                            lhs.sequence,
                            reference[lhs.end:end]])
    rhs_observed = "".join([reference[start:rhs.start],
                            rhs.sequence,
                            reference[rhs.end:end]])
    return sequence_based_contains(reference[start:end], lhs_observed, rhs_observed)


def is_contained(reference, lhs, rhs):
    """Check if `lhs` is contained in `rhs`."""
    return contains(reference, rhs, lhs)


def are_disjoint(reference, lhs, rhs):
    """Check if two variants are disjoint."""
    if lhs == rhs:
        return False
    if not lhs or not rhs or lhs.is_disjoint(rhs):
        return True

    return graph_based_are_disjoint(reference, LCSgraph.from_supremal(reference, lhs), LCSgraph.from_supremal(reference, rhs))


def have_overlap(reference, lhs, rhs):
    """Check if two variants overlap."""
    if lhs == rhs or not lhs or not rhs or lhs.is_disjoint(rhs):
        return False

    return graph_based_have_overlap(reference, LCSgraph.from_supremal(reference, lhs), LCSgraph.from_supremal(reference, rhs))


def compare(reference, lhs, rhs):
    """Compare two supremal variants.

    Parameters
    ----------
    reference : str
        The reference sequence.
    lhs : `Variant`
        The supremal variant on the left-hand side.
    rhs : `Variant`
        The supremal variant on the right-hand side.

    Returns
    -------
    `Relation`
        The relation between the two supremal variants.
    """

    if lhs == rhs:
        return Relation.EQUIVALENT

    if not lhs or not rhs or lhs.is_disjoint(rhs):
        return Relation.DISJOINT

    return graph_based_compare(reference, LCSgraph.from_supremal(reference, lhs), LCSgraph.from_supremal(reference, rhs))
