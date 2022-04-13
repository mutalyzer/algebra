import sys
from algebra.lcs.onp import edit as edit_onp
from algebra.lcs.wupp import edit as edit_test, lcs_graph as graph_test, traversal
from algebra.lcs.efficient import edit as edit_gold, build as graph_gold, traversal as traversal_gold
from algebra.variants.variant import to_hgvs, Variant, patch, merge_co_insertions, turbo_sort
from pprint import pprint
import random


def compare_matrix(test, gold, f):
    assert len(test) == len(gold)
    assert len(test[0]) == len(gold[0])
    row_len = len(test)
    col_len = len(test[0])

    for row in range(row_len):
        for col in range(col_len):
            if gold[row][col] is not None and gold[row][col] <= f:
                assert gold[row][col] == test[row][col]


def main():
    min_rand = 10
    max_rand = 10

    if len(sys.argv) == 1:
        reference = "".join(random.choice("ACGT") for _ in range(random.randint(min_rand, max_rand)))
        observed = "".join(random.choice("ACGT") for _ in range(random.randint(min_rand, max_rand)))
    elif len(sys.argv) == 3:
        reference = sys.argv[1]
        observed = sys.argv[2]
    else:
        raise("usage")

    print(reference, observed)

    dist_gold, nodes_gold = edit_gold(reference, observed)
    print(dist_gold)
    dist_test, matrix_test, nodes_test = edit_test(reference, observed)
    print(dist_test)
    # pprint(matrix_test)
    # pprint(matrix_gold)

    assert dist_test == dist_gold

    # compare_matrix(matrix_test, matrix_gold, dist_test)

    print("gold")
    _, lcs_graph_gold = graph_gold(nodes_gold, reference, observed)
    paths_gold = traversal_gold(reference, observed, lcs_graph_gold, atomics=True)
    hgvs_gold = []
    for path in paths_gold:
        hgvs_gold.append(to_hgvs(path, reference))
    print(f"length {len(hgvs_gold)}/{len(set(hgvs_gold))}")
    for h in sorted(hgvs_gold):
        print(h)

    print("test")
    #print(nodes_test)
    for level in nodes_test:
        print(level)

    graph = graph_test(reference, observed, nodes_test)
    print("digraph {")
    for node, edge_list in graph.items():
        print(node)
        for child, edge in edge_list:
            # print(node, child, to_hgvs(edge, reference))
            # us "0_0" -> "1_1" [label="="];
            print(f'"{node[0]}_{node[1]}" -> "{child[0]}_{child[1]}" [label="{to_hgvs(edge, reference)}"];')
    print("}")

    paths_test = traversal(reference, observed, graph, atomics=True)
    # s = set()
    hgvs_test = []
    for path in paths_test:
        var = merge_co_insertions(turbo_sort(path))
        # for v in path:
        #     print(v.to_hgvs(reference))
        # print(to_hgvs(path, reference, sort=False))
        # print(to_hgvs(var, reference, sort=False))
        # print()
        hgvs_test.append(to_hgvs(var, reference, sort=False))
        # s.add(hgvs)
    print(f"length {len(hgvs_test)}/{len(set(hgvs_test))}")
    for h in sorted(hgvs_test):
        print(h)
    # assert len(s) == len(paths)

    # assert set(hgvs_test) == set(hgvs_gold)
    print(set(hgvs_test) - set(hgvs_gold))


if __name__ == '__main__':
    main()
