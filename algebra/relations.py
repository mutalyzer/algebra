from enum import Enum
from itertools import product
from .lcs.onp import edit as edit_distance_only
from .lcs.wupp import edit, lcs_graph


class Relation(Enum):
    EQUIVALENT = "equivalent"
    CONTAINS = "contains"
    IS_CONTAINED = "is_contained"
    OVERLAP = "overlap"
    DISJOINT = "disjoint"


def are_equivalent(_reference, lhs, rhs):
    return lhs == rhs


def contains(reference, lhs, rhs):
    if lhs == rhs:
        return False

    lhs_distance = edit_distance_only(reference, lhs)
    rhs_distance = edit_distance_only(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if lhs_distance + rhs_distance == distance:
        return False

    return lhs_distance - rhs_distance == distance


def is_contained(reference, lhs, rhs):
    if lhs == rhs:
        return False

    lhs_distance = edit_distance_only(reference, lhs)
    rhs_distance = edit_distance_only(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if lhs_distance + rhs_distance == distance:
        return False

    return rhs_distance - lhs_distance == distance


def are_disjoint(reference, lhs, rhs):
    if lhs == rhs:
        return False

    lhs_distance, lhs_lcs_nodes = edit(reference, lhs)
    rhs_distance, rhs_lcs_nodes = edit(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if lhs_distance + rhs_distance == distance:
        return True

    if (lhs_distance - rhs_distance == distance or
            rhs_distance - lhs_distance == distance):
        return False

    _, lhs_edges = lcs_graph(reference, lhs, lhs_lcs_nodes)
    _, rhs_edges = lcs_graph(reference, rhs, rhs_lcs_nodes)

    for lhs_variant, rhs_variant in product(lhs_edges, rhs_edges):
        if not lhs_variant.is_disjoint(rhs_variant):
            return False

    return True


def have_overlap(reference, lhs, rhs):
    if lhs == rhs:
        return False

    lhs_distance, lhs_lcs_nodes = edit(reference, lhs)
    rhs_distance, rhs_lcs_nodes = edit(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if (lhs_distance + rhs_distance == distance or
            lhs_distance - rhs_distance == distance or
            rhs_distance - lhs_distance == distance):
        return False

    _, lhs_edges = lcs_graph(reference, lhs, lhs_lcs_nodes)
    _, rhs_edges = lcs_graph(reference, rhs, rhs_lcs_nodes)

    for lhs_variant, rhs_variant in product(lhs_edges, rhs_edges):
        if not lhs_variant.is_disjoint(rhs_variant):
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

    for lhs_variant, rhs_variant in product(lhs_edges, rhs_edges):
        if not lhs_variant.is_disjoint(rhs_variant):
            return Relation.OVERLAP

    return Relation.DISJOINT
