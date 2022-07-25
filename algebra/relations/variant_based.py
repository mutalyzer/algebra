from .supremal_based import compare as compare_supremal, find_supremal, spanning_variant
from ..variants import patch


def compare(reference, lhs, rhs):
    lhs_sup = find_supremal(reference, spanning_variant(reference, patch(reference, lhs), lhs))
    rhs_sup = find_supremal(reference, spanning_variant(reference, patch(reference, rhs), rhs))

    return compare_supremal(reference, lhs_sup, rhs_sup)

