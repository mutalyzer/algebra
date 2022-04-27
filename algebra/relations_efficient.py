from .lcs.efficient import edit, build
from .lcs.onp import edit as edit_distance_only
from .relations import Relation


def are_disjoint(reference, lhs, rhs):
    if lhs == rhs:
        return False

    lhs_distance, lhs_lcs = edit(reference, lhs)
    rhs_distance, rhs_lcs = edit(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if (lhs_distance - rhs_distance == distance or
            rhs_distance - lhs_distance == distance):
        return False

    if lhs_distance + rhs_distance == distance:
        return True

    lhs_ops, _ = build(lhs_lcs, reference, lhs)
    rhs_ops, _ = build(rhs_lcs, reference, rhs)

    if lhs_ops.isdisjoint(rhs_ops):
        return True

    return False


def have_overlap(reference, lhs, rhs):
    if lhs == rhs:
        return False

    lhs_distance, lhs_lcs = edit(reference, lhs)
    rhs_distance, rhs_lcs = edit(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if (lhs_distance - rhs_distance == distance or
            rhs_distance - lhs_distance == distance or
            lhs_distance + rhs_distance == distance):
        return False

    lhs_ops, _ = build(lhs_lcs, reference, lhs)
    rhs_ops, _ = build(rhs_lcs, reference, rhs)

    if lhs_ops.isdisjoint(rhs_ops):
        return True

    return False


def compare(reference, lhs, rhs):
    if lhs == rhs:
        return Relation.EQUIVALENT

    lhs_distance, lhs_lcs = edit(reference, lhs)
    rhs_distance, rhs_lcs = edit(reference, rhs)
    distance = edit_distance_only(lhs, rhs)

    if lhs_distance + rhs_distance == distance:
        return Relation.DISJOINT

    if lhs_distance - rhs_distance == distance:
        return Relation.CONTAINS

    if rhs_distance - lhs_distance == distance:
        return Relation.IS_CONTAINED

    lhs_ops, _ = build(lhs_lcs, reference, lhs)
    rhs_ops, _ = build(rhs_lcs, reference, rhs)

    if lhs_ops.isdisjoint(rhs_ops):
        return Relation.DISJOINT

    return Relation.OVERLAP
