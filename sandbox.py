from itertools import product
import sys
from algebra.lcs.wupp import edit, lcs_graph, to_dot
from algebra.utils import random_sequence
from algebra.variants.variant import Variant


def are_disjoint(lhs, rhs):
    for lhs_variant, rhs_variant in product(lhs, rhs):
        if not lhs_variant.is_disjoint(rhs_variant):
            return False
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


def compare(reference, lhs, rhs):
    print(f'"{reference}" "{lhs}" "{rhs}"')

    lhs_distance, lhs_lcs_nodes = edit(reference, lhs)
    rhs_distance, rhs_lcs_nodes = edit(reference, rhs)

    lhs_root, lhs_edges = lcs_graph(reference, lhs, lhs_lcs_nodes)
    rhs_root, rhs_edges = lcs_graph(reference, rhs, rhs_lcs_nodes)

    print(to_dot(reference, lhs_root))
    print(to_dot(reference, rhs_root))

    lhs_ops = ops_set(lhs_edges)
    rhs_ops = ops_set(rhs_edges)

    print(sorted({variant.to_hgvs(reference) for variant in lhs_ops}))
    print(sorted({variant.to_hgvs(reference) for variant in rhs_ops}))

    print("EDGES")
    print([variant.to_hgvs(reference) for variant in lhs_edges])
    print([variant.to_hgvs(reference) for variant in rhs_edges])
    print("-----")

    print(lhs_ops.isdisjoint(rhs_ops))
    print(are_disjoint(lhs_edges, rhs_edges))

    assert lhs_ops.isdisjoint(rhs_ops) == are_disjoint(lhs_edges, rhs_edges)


def main():
    if len(sys.argv) < 4:
        reference = random_sequence(15)
        lhs = random_sequence(15)
        rhs = random_sequence(15)
    else:
        reference = sys.argv[1]
        lhs = sys.argv[2]
        rhs = sys.argv[3]

    compare(reference, lhs, rhs)


if __name__ == "__main__":
    main()
