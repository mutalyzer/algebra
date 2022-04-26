from itertools import product
import sys
from algebra.lcs.wupp import edit, lcs_graph, traversal, to_dot
from algebra.relations import disjoint_variants, ops_set
from algebra.utils import random_sequence


def are_disjoint(lhs, rhs):
    for lhs_variant, rhs_variant in product(lhs, rhs):
        if not disjoint_variants(lhs_variant, rhs_variant):
            print(lhs_variant.to_hgvs(), rhs_variant.to_hgvs())
            return False
    return True


def main():
    if len(sys.argv) < 4:
        reference = random_sequence(15)
        lhs = random_sequence(15)
        rhs = random_sequence(15)
    else:
        reference = sys.argv[1]
        lhs = sys.argv[2]
        rhs = sys.argv[3]

    print(f'"{reference}" "{lhs}" "{rhs}"')

    lhs_distance, lhs_lcs_nodes = edit(reference, lhs)
    rhs_distance, rhs_lcs_nodes = edit(reference, rhs)

    lhs_root, lhs_edges = lcs_graph(reference, lhs, lhs_lcs_nodes)
    rhs_root, rhs_edges = lcs_graph(reference, rhs, rhs_lcs_nodes)

    print(to_dot(reference, lhs_root))
    print(to_dot(reference, rhs_root))

    print(sorted([variant.to_hgvs(reference) for variant in lhs_edges]))
    print(sorted([variant.to_hgvs(reference) for variant in rhs_edges]))

    lhs_ops = ops_set(lhs_edges)
    rhs_ops = ops_set(rhs_edges)

    print(sorted({variant.to_hgvs(reference) for variant in lhs_ops}))
    print(sorted({variant.to_hgvs(reference) for variant in rhs_ops}))

    print(lhs_ops.isdisjoint(rhs_ops))
    print(are_disjoint(lhs_edges, rhs_edges))


if __name__ == "__main__":
    main()
