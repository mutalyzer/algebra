from .lcs.efficient import edit, lcs_graph
from .lcs.onp import edit as edit_fast


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

    # FIXME: disjoint check
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

    # FIXME: disjoint check
    return True


def compare(reference, lhs, rhs):
    if lhs == rhs:
        return "equivalent"

    lhs_distance, lhs_lcs = edit(reference, lhs)
    rhs_distance, rhs_lcs = edit(reference, rhs)
    distance = edit_fast(lhs, rhs)

    if lhs_distance + rhs_distance == distance:
        return "disjoint"

    if lhs_distance - rhs_distance == distance:
        return "contains"

    if rhs_distance - lhs_distance == distance:
        return "is_contained"

    lhs_graph = lcs_graph(reference, lhs, lhs_lcs)
    rhs_graph = lcs_graph(reference, rhs, rhs_lcs)

    # FIXME: disjoint check
    return "overlap"
