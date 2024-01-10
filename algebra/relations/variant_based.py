"""Functions to compare variant alleles."""


from ..lcs import LCSgraph
from .relation import Relation
from .graph_based import compare as compare_graph


def are_equivalent(reference, lhs, rhs):
    """Check if two variants are equivalent."""
    return compare(reference, lhs, rhs) == Relation.EQUIVALENT


def contains(reference, lhs, rhs):
    """Check if `lhs` contains `rhs`."""
    return compare(reference, lhs, rhs) == Relation.CONTAINS


def is_contained(reference, lhs, rhs):
    """Check if `lhs` is contained in `rhs`."""
    return compare(reference, lhs, rhs) == Relation.IS_CONTAINED


def are_disjoint(reference, lhs, rhs):
    """Check if two variants are disjoint."""
    return compare(reference, lhs, rhs) == Relation.DISJOINT


def have_overlap(reference, lhs, rhs):
    """Check if two variants overlap."""
    return compare(reference, lhs, rhs) == Relation.OVERLAP


def compare(reference, lhs, rhs):
    """Compare two variant alleles.

    Parameters
    ----------
    reference : str
        The reference sequence.
    lhs : list of variants
        The variant allele on the left-hand side.
    rhs : list of variants
        The variant allele on the right-hand side.

    Returns
    -------
    `Relation`
        The relation between the two variant alleles.
    """

    return compare_graph(reference, LCSgraph.from_variant(reference, lhs), LCSgraph.from_variant(reference, rhs))
