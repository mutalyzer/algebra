"""Functions to find supremal LCS graphs."""


from operator import attrgetter
from os.path import commonprefix
from ..lcs.lcs_graph import LCSgraph
from ..variants import Variant, patch


def lcs_graph(reference, variants, offset=10):
    """Iteratively find the supremal LCS graph for an allele by repeatedly
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
    graph : `LCSgraph`
        The LCS graph.
    """

    if not variants:
        return LCSgraph("", "")

    start = min(variants, key=attrgetter("start")).start
    end = max(variants, key=attrgetter("end")).end
    observed = patch(reference[start:end], [Variant(variant.start - start, variant.end - start, variant.sequence) for variant in variants])
    variant = Variant(start, end, observed)

    offset = max(offset, len(variant) // 2, 1)

    while True:
        start = max(0, variant.start - offset)
        end = min(len(reference), variant.end + offset)
        observed = reference[start:variant.start] + variant.sequence + reference[variant.end:end]

        graph = LCSgraph(reference[start:end], observed, shift=start)

        if not graph.supremal:
            return graph

        if ((graph.supremal.start > start or graph.supremal.start == 0) and
                (graph.supremal.end < end or graph.supremal.end == len(reference))):
            return graph

        offset *= 2


def lcs_graph_sequence(reference, observed):
    """The supremal LCS graph for two sequences."""
    if reference == observed:
        return LCSgraph("", "")

    prefix_len, suffix_len = trim(reference, observed)
    return lcs_graph(reference, [Variant(prefix_len, len(reference) - suffix_len, observed[prefix_len:len(observed) - suffix_len])])


def lcs_graph_supremal(reference, supremal):
    """The supremal LCS graph for a supremal variant."""
    return LCSgraph(reference[supremal.start:supremal.end], supremal.sequence, shift=supremal.start)


def trim(lhs, rhs):
    """Find the lengths of the common prefix and common suffix between
    two sequences."""
    idx = len(commonprefix([lhs, rhs]))
    return idx, len(commonprefix([lhs[idx:][::-1], rhs[idx:][::-1]]))
