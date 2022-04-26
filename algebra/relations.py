from enum import Enum
from itertools import product
from .lcs.onp import edit as edit_distance_only
from .lcs.wupp import edit, lcs_graph
from .variants.variant import Variant


class Relation(Enum):
    EQUIVALENT = "equivalent"
    CONTAINS = "contains"
    IS_CONTAINED = "is_contained"
    OVERLAP = "overlap"
    DISJOINT = "disjoint"


def disjoint_variants(lhs, rhs):
    if lhs.start < rhs.end and rhs.start < lhs.end and lhs.start < lhs.end and rhs.start < rhs.end:
        return False
    if lhs.start <= rhs.end and rhs.start <= lhs.end:
        return set(lhs.sequence).isdisjoint(set(rhs.sequence))
    return True


def ops_set(edges):
    def explode(variant):
        for pos in range(variant.start, variant.end):
            yield Variant(pos, pos + 1)
        for pos in range(variant.start, variant.end + 1):
            for symbol in variant.sequence:
                yield Variant(pos, pos, symbol)

    ops = set()
    for edge in edges:
        ops.update(explode(edge))

    return ops


def are_equivalent(reference, lhs, rhs):
    return lhs == rhs


def contains(reference, lhs, rhs):
    return (lhs != rhs and
            edit_distance_only(reference, lhs) - edit_distance_only(reference, rhs) == edit_distance_only(lhs, rhs))


def is_contained(reference, lhs, rhs):
    return (lhs != rhs and
            edit_distance_only(reference, rhs) - edit_distance_only(reference, lhs) == edit_distance_only(lhs, rhs))


def are_disjoint(reference, lhs, rhs):
    lhs_distance, lhs_lcs_nodes = edit(reference, lhs)
    rhs_distance, rhs_lcs_nodes = edit(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if (lhs == rhs or
            lhs_distance - rhs_distance == distance or
            rhs_distance - lhs_distance == distance):
        return False

    if lhs_distance + rhs_distance == distance:
        return True

    _, lhs_edges = lcs_graph(reference, lhs, lhs_lcs_nodes)
    _, rhs_edges = lcs_graph(reference, rhs, rhs_lcs_nodes)

    lhs_ops = ops_set(lhs_edges)
    rhs_ops = ops_set(rhs_edges)

    if lhs_ops.isdisjoint(rhs_ops):
        return True

    return False


def have_overlap(reference, lhs, rhs):
    lhs_distance, lhs_lcs_nodes = edit(reference, lhs)
    rhs_distance, rhs_lcs_nodes = edit(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if (lhs == rhs or
            lhs_distance - rhs_distance == distance or
            rhs_distance - lhs_distance == distance or
            lhs_distance + rhs_distance == distance):
        return False

    _, lhs_edges = lcs_graph(reference, lhs, lhs_lcs_nodes)
    _, rhs_edges = lcs_graph(reference, rhs, rhs_lcs_nodes)

    lhs_ops = ops_set(lhs_edges)
    rhs_ops = ops_set(rhs_edges)

    if lhs_ops.isdisjoint(rhs_ops):
        return True

    return False


def compare(reference, lhs, rhs):
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

    # lhs_ops = ops_set(lhs_edges)
    # rhs_ops = ops_set(rhs_edges)
    #
    # if lhs_ops.isdisjoint(rhs_ops):
    #     return Relation.DISJOINT
    #
    # return Relation.OVERLAP

    for lhs_edge, rhs_edge in product(lhs_edges, rhs_edges):
        if not disjoint_variants(lhs_edge, rhs_edge):
            return Relation.OVERLAP

    return Relation.DISJOINT
