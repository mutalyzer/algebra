import sys
from algebra.lcs.efficient import edit as edit_old, build, traversal as traversal_old
from algebra.lcs.wupp import edit as edit_new, lcs_graph, traversal as traversal_new
from algebra.relations import disjoint_variants
from algebra.variants.variant import Variant, to_hgvs


def main():
    if len(sys.argv) < 3:
        reference = "CATATATCG"
        observed = "CTTATAGCAT"
    else:
        reference = sys.argv[1]
        observed = sys.argv[2]

    distance_old, lcs_nodes_old = edit_old(reference, observed)
    ops_old, graph_old = build(lcs_nodes_old, reference, observed)

    hgvs = set()
    for op in ops_old:
        if op[1] == "ins":
            variant = Variant(op[0], op[0], op[2])
        else:
            variant = Variant(op[0] - 1, op[0])
        hgvs.add(variant.to_hgvs())

    print(sorted(hgvs))
    print(len(hgvs))

    distance_new, lcs_nodes_new = edit_new(reference, observed)
    #ops_new = ops_set(reference, observed, lcs_nodes_new)
    graph_new, edges = lcs_graph(reference, observed, lcs_nodes_new)

    #print(ops_new)
    print(edges)
    for i in range(len(edges)):
        for j in range(i, len(edges)):
            print(edges[i].to_hgvs(reference), edges[j].to_hgvs(reference), disjoint_variants(edges[i], edges[j]))


if __name__ == "__main__":
    main()
