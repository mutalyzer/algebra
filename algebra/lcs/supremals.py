"""Functions to find supremal variants."""


from operator import attrgetter
from os.path import commonprefix
from ..lcs.all_lcs import LCSnode, bfs_traversal, lcs_graph
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
    source : `LCS_Node`
        The source of the LCS graph in which the supremal was determined.
    """

    if not variants:
        return Variant(0, 0, ""), LCSnode(0, 0, len(reference))

    start = min(variants, key=attrgetter("start")).start
    end = max(variants, key=attrgetter("end")).end
    observed = patch(reference[start:end], [Variant(variant.start - start, variant.end - start, variant.sequence) for variant in variants])
    variant = Variant(start, end, observed)

    offset = max(offset, len(variant) // 2, 1)

    while True:
        start = max(0, variant.start - offset)
        end = min(len(reference), variant.end + offset)
        observed = reference[start:variant.start] + variant.sequence + reference[variant.end:end]

        graph = lcs_graph(reference[start:end], observed, shift=start)
        edges = {edge[0] for *_, edge in bfs_traversal(graph)}

        if not edges:
            return Variant(0, 0, ""), graph

        edges_start = min(edges, key=attrgetter("start")).start
        edges_end = max(edges, key=attrgetter("end")).end

        if (edges_start > start or edges_start == 0) and (edges_end < end or edges_end == len(reference)):
            return Variant(edges_start, edges_end, observed[edges_start - start:len(observed) - (end - edges_end)]), graph

        offset *= 2


def supremal_sequence(reference, observed):
    """The supremal variant for two sequences."""
    if reference == observed:
        return supremal(reference, [])

    prefix_len, suffix_len = trim(reference, observed)
    return supremal(reference, [Variant(prefix_len, len(reference) - suffix_len, observed[prefix_len:len(observed) - suffix_len])])
