import sys
from algebra.lcs.wupp import edit, lcs_graph, traversal
from algebra.relations import disjoint_variants, ops_set
from algebra.utils import random_sequence


def main():
    if len(sys.argv) < 4:
        reference = random_sequence(15)
        lhs = random_sequence(15)
        rhs = random_sequence(15)
    else:
        reference = sys.argv[1]
        lhs = sys.argv[2]
        rhs = sys.argv[3]

    print(reference, lhs, rhs)

    lhs_distance, lhs_lcs_nodes = edit(reference, lhs)
    rhs_distance, rhs_lcs_nodes = edit(reference, rhs)

    lhs_root, lhs_edges = lcs_graph(reference, lhs, lhs_lcs_nodes)
    rhs_root, rhs_edges = lcs_graph(reference, rhs, rhs_lcs_nodes)

    print({var.to_hgvs(reference) for var in ops_set(lhs_edges)})
    print({var.to_hgvs(reference) for var in ops_set(rhs_edges)})



if __name__ == "__main__":
    main()
