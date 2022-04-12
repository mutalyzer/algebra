from .lcs.efficient import edit, build
from .lcs.onp import edit as edit_fast
from enum import Enum


class Relation(Enum):
    EQUIVALENT = "equivalent"
    CONTAINS = "contains"
    IS_CONTAINED = "is_contained"
    OVERLAP = "overlap"
    DISJOINT = "disjoint"


def are_equivalent(reference, lhs, rhs):
    return lhs == rhs


def are_disjoint(reference, lhs, rhs):
    lhs_distance, lhs_lcs = edit(reference, lhs)
    rhs_distance, rhs_lcs = edit(reference, rhs)
    distance = edit_fast(lhs, rhs)

    if (lhs == rhs or
            lhs_distance - rhs_distance == distance or
            rhs_distance - lhs_distance == distance):
        return False

    if lhs_distance + rhs_distance == distance:
        return True

    lhs_ops, lhs_graph = build(lhs_lcs, reference, lhs)
    rhs_ops, rhs_graph = build(rhs_lcs, reference, rhs)

    if lhs_ops.isdisjoint(rhs_ops):
        return True

    return False


def contains(reference, lhs, rhs):
    return (lhs != rhs and
            edit_fast(reference, lhs) - edit_fast(reference, rhs) == edit_fast(lhs, rhs))


def is_contained(reference, lhs, rhs):
    return (lhs != rhs and
            edit_fast(reference, rhs) - edit_fast(reference, lhs) == edit_fast(lhs, rhs))


def have_overlap(reference, lhs, rhs):
    lhs_distance, lhs_lcs = edit(reference, lhs)
    rhs_distance, rhs_lcs = edit(reference, rhs)
    distance = edit_fast(lhs, rhs)

    if (lhs == rhs or
            lhs_distance - rhs_distance == distance or
            rhs_distance - lhs_distance == distance or
            lhs_distance + rhs_distance == distance):
        return False

    lhs_ops, lhs_graph = build(lhs_lcs, reference, lhs)
    rhs_ops, rhs_graph = build(rhs_lcs, reference, rhs)

    if lhs_ops.isdisjoint(rhs_ops):
        return True

    return False


def compare(reference, lhs, rhs):
    if lhs == rhs:
        return Relation.EQUIVALENT

    lhs_distance, lhs_lcs = edit(reference, lhs)
    rhs_distance, rhs_lcs = edit(reference, rhs)
    distance = edit_fast(lhs, rhs)

    if lhs_distance + rhs_distance == distance:
        return Relation.DISJOINT

    if lhs_distance - rhs_distance == distance:
        return Relation.CONTAINS

    if rhs_distance - lhs_distance == distance:
        return Relation.IS_CONTAINED

    lhs_ops, lhs_graph = build(lhs_lcs, reference, lhs)
    rhs_ops, rhs_graph = build(rhs_lcs, reference, rhs)

    if lhs_ops.isdisjoint(rhs_ops):
        return Relation.DISJOINT

    return Relation.OVERLAP
