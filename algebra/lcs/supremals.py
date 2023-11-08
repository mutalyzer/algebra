"""Functions to find supremal variants."""


from operator import attrgetter
from os.path import commonprefix
from ..lcs.all_lcs import _Node, build_graph, edit
from ..variants import Variant, patch


def trim(lhs, rhs):
    """Find the lengths of the common prefix and common suffix between
    two sequences."""
    idx = len(commonprefix([lhs, rhs]))
    return idx, len(commonprefix([lhs[idx:][::-1], rhs[idx:][::-1]]))


def supremal(reference, variants, offset=10):
    """Iteratively find the supremal variant for an allele by repeatedly
    widening a range of influence.

    Parameters
    ----------
    reference : str
        The reference sequence.
    variants : list of `Variant`s
        The allele of interest.

    Other Parameters
    ----------------
    offset : int, optional
        The minimum offset around the variant.

    Returns
    -------
    supremal : `Variant`
        The supremal variant.
    root : `_Node`
        The LCS graph in which the supremal was determined.
    observed : str
        The observed sequence for the given LCS graph.
    shift : int
        The offset by which the variants in the LCS graph are shifted.
    """

    def spanning_variant(reference, observed, variants, shift=0):
        if not variants:
            return Variant(0, 0, "")

        start = min(variants, key=attrgetter("start")).start
        end = max(variants, key=attrgetter("end")).end
        return Variant(start, end, observed[shift + start:shift + len(observed) - (len(reference) - end)])

    observed = patch(reference, variants)
    variant = spanning_variant(reference, observed, variants)
    if not variant:
        return variant, _Node(0, 0, len(reference)), observed, 0

    offset = max(offset, len(variant) // 2, 1)

    while True:
        start = max(0, variant.start - offset)
        end = min(len(reference), variant.end + offset)

        observed = reference[start:variant.start] + variant.sequence + reference[variant.end:end]

        _, lcs_nodes = edit(reference[start:end], observed)
        root, edges = build_graph(reference[start:end], observed, lcs_nodes, start)
        variant = spanning_variant(reference[start:end], observed, edges, -start)

        if not variant:
            return variant, root, observed, start

        if ((variant.start > start or variant.start == 0) and
                (variant.end < end or variant.end == len(reference))):
            return variant, root, observed, start

        offset *= 2


def supremal_sequence(reference, observed):
    """The supremal variant for two sequences."""
    if reference == observed:
        return supremal(reference, [])

    prefix_len, suffix_len = trim(reference, observed)
    return supremal(reference, [Variant(prefix_len, len(reference) - suffix_len, observed[prefix_len:len(observed) - suffix_len])])
