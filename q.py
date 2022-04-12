import sys
from algebra.lcs.onp import edit as edit_onp
from algebra.lcs.wupp import edit as edit_test, lcs_graph as graph_test, traversal
from algebra.lcs.efficient import edit as edit_gold, build as graph_gold, traversal as traversal_gold
from algebra.variants.variant import to_hgvs, Variant, patch, merge_cons
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
    min_rand = 1
    max_rand = 5

    if len(sys.argv) == 1:
        reference = "".join(random.choice("ACGT") for _ in range(random.randint(min_rand, max_rand)))
        observed = "".join(random.choice("ACGT") for _ in range(random.randint(min_rand, max_rand)))
    elif len(sys.argv) == 3:
        reference = sys.argv[1]
        observed = sys.argv[2]
    else:
        raise("usage")

    print(reference, observed)

    dist_test, matrix_test, nodes_test = edit_test(reference, observed)
    dist_gold, nodes_gold = edit_gold(reference, observed)
    print(dist_test)
    # pprint(matrix_test)
    print(dist_gold)
    # pprint(matrix_gold)

    assert dist_test == dist_gold

    # compare_matrix(matrix_test, matrix_gold, dist_test)

    print(nodes_test)
    for level in nodes_test:
        print(level)

    graph = graph_test(reference, observed, nodes_test)

    print("test")
    paths_test = traversal(reference, observed, graph, atomics=True)
    # s = set()
    hgvs_test = []
    for path in paths_test:
        hgvs_test.append(to_hgvs(path, reference))
        # s.add(hgvs)
    for h in sorted(hgvs_test):
        print(h)
    # assert len(s) == len(paths)

    # reference = 'ABBC'
    # v = Variant(1, 2, 'X')
    # print(v.to_hgvs(reference))
    # for atomic in v.atomics():
    #     print(to_hgvs(atomic, reference), patch(reference, atomic))

    print("gold")
    _, lcs_graph_gold = graph_gold(nodes_gold, reference, observed)
    paths_gold = traversal_gold(reference, observed, lcs_graph_gold, atomics=True)
    hgvs_gold = []
    for path in paths_gold:
        hgvs_gold.append(to_hgvs(path, reference))
    for h in sorted(hgvs_gold):
        print(h)

    assert set(hgvs_test) == set(hgvs_gold)


if __name__ == '__main__':
    main()
