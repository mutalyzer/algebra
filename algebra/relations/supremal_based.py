"""Functions to find and compare supremal (minimum spanning) variants."""


from operator import attrgetter
from .relation import Relation
from .sequence_based import (are_disjoint as sequence_based_are_disjoint,
                             compare as sequence_based_compare,
                             contains as sequence_based_contains,
                             have_overlap as sequence_based_have_overlap)
from ..lcs import edit, lcs_graph
from ..variants import Variant


def are_equivalent(_reference, lhs, rhs):
    """Check if two variants are equivalent."""
    return lhs == rhs


def contains(reference, lhs, rhs):
    """Check if `lhs` contains `rhs`."""
    start = min(lhs.start, rhs.start)
    end = max(lhs.end, rhs.end)
    lhs_observed = reference[min(start, lhs.start):lhs.start] + lhs.sequence + reference[lhs.end:max(end, lhs.end)]
    rhs_observed = reference[min(start, rhs.start):rhs.start] + rhs.sequence + reference[rhs.end:max(end, rhs.end)]
    return sequence_based_contains(reference[start:end], lhs_observed, rhs_observed)


def is_contained(reference, lhs, rhs):
    """Check if `lhs` is contained in `rhs`."""
    return contains(reference, rhs, lhs)


def are_disjoint(reference, lhs, rhs):
    """Check if two variants are disjoint."""
    if lhs.is_disjoint(rhs):
        return True

    start = min(lhs.start, rhs.start)
    end = max(lhs.end, rhs.end)
    lhs_observed = reference[min(start, lhs.start):lhs.start] + lhs.sequence + reference[lhs.end:max(end, lhs.end)]
    rhs_observed = reference[min(start, rhs.start):rhs.start] + rhs.sequence + reference[rhs.end:max(end, rhs.end)]
    return sequence_based_are_disjoint(reference[start:end], lhs_observed, rhs_observed)


def have_overlap(reference, lhs, rhs):
    """Check if two variants overlap."""
    start = min(lhs.start, rhs.start)
    end = max(lhs.end, rhs.end)
    lhs_observed = reference[min(start, lhs.start):lhs.start] + lhs.sequence + reference[lhs.end:max(end, lhs.end)]
    rhs_observed = reference[min(start, rhs.start):rhs.start] + rhs.sequence + reference[rhs.end:max(end, rhs.end)]
    return sequence_based_have_overlap(reference[start:end], lhs_observed, rhs_observed)


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

    if lhs.is_disjoint(rhs):
        return Relation.DISJOINT

    start = min(lhs.start, rhs.start)
    end = max(lhs.end, rhs.end)
    lhs_observed = reference[min(start, lhs.start):lhs.start] + lhs.sequence + reference[lhs.end:max(end, lhs.end)]
    rhs_observed = reference[min(start, rhs.start):rhs.start] + rhs.sequence + reference[rhs.end:max(end, rhs.end)]
    return sequence_based_compare(reference[start:end], lhs_observed, rhs_observed)


def spanning_variant(reference, observed, variants):
    """Calculate the minimum spanning variant for a collection of
    variants. If the collection of variants is the collection of all
    minimal variants the minimum spanning variant is the supremal
    variant.

    See Also
    --------
    `algebra.lcs.all_lcs.lcs_graph` : The collection of all minimal variants.
    """

    if not variants:
        return Variant(0, 0, "")

    start = min(variants, key=attrgetter("start")).start
    end = max(variants, key=attrgetter("end")).end
    return Variant(start, end, observed[start:len(observed) - (len(reference) - end)])


def find_supremal(reference, variant, offset=10):
    """Iteratively find the supremal variant.

    Parameters
    ----------
    reference : str
        The reference sequence.
    variant : `Variant`
        The variant of interest. Allele descriptions should be converted
        to a single (delins type) variant.

    Other Parameters
    ----------------
    offset : int, optional
        The minimum offset around the variant.

    See Also
    --------
    `spanning_variant` : The minimum spanning (delins) variant.
    """

    offset = max(offset, len(variant), 1)

    while True:
        start = max(0, variant.start - offset)
        end = min(len(reference), variant.end + offset)

        observed = reference[start:variant.start] + variant.sequence + reference[variant.end:end]

        _, lcs_nodes = edit(reference[start:end], observed)
        _, edges = lcs_graph(reference[start:end], observed, lcs_nodes)
        supremum = spanning_variant(reference[start:end], observed, edges)

        supremum.start += start
        supremum.end += start

        if ((supremum.start > start or supremum.start == 0) and
                (supremum.end < end or supremum.end == len(reference))):
            return supremum

        offset *= 2
