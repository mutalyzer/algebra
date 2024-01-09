from itertools import product
from ..lcs import edit_distance
from .relation import Relation


def compare(reference, lhs, rhs):
    if lhs.supremal == rhs.supremal:
        return Relation.EQUIVALENT

    if lhs.supremal.is_disjoint(rhs.supremal):
        return Relation.DISJOINT

    start = min(lhs.supremal.start, rhs.supremal.start)
    end = max(lhs.supremal.end, rhs.supremal.end)

    # TODO: class method on LCSgraph?
    def observed(supremal):
        return (reference[min(start, supremal.start):supremal.start] +
                supremal.sequence +
                reference[supremal.end:max(end, supremal.end)])

    lhs_observed = observed(lhs.supremal)
    rhs_observed = observed(rhs.supremal)
    distance = edit_distance(lhs_observed, rhs_observed)

    if lhs.distance + rhs.distance == distance:
        return Relation.DISJOINT

    if lhs.distance - rhs.distance == distance:
        return Relation.CONTAINS

    if rhs.distance - lhs.distance == distance:
        return Relation.IS_CONTAINED

    for lhs_variant, rhs_variant in product(lhs.edges(), rhs.edges()):
        if not lhs_variant.is_disjoint(rhs_variant):
            return Relation.OVERLAP

    return Relation.DISJOINT
