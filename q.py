import sys
from algebra.lcs.onp import edit as edit_onp
from algebra.lcs.wupp import edit as edit_test, lcs_graph as graph_test, traversal, to_dot as to_dot_test
from algebra.lcs.efficient import edit as edit_gold, build as graph_gold, traversal as traversal_gold, to_dot as to_dot_gold
from algebra.variants.variant import to_hgvs, Variant, patch
from pprint import pprint
import random


def main():
    min_rand = 1
    max_rand = 15

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
    dist_test, nodes_test = edit_test(reference, observed)
    print(dist_test)

    assert dist_test == dist_gold

    print("gold")
    _, lcs_graph_gold = graph_gold(nodes_gold, reference, observed)
    print(to_dot_gold(reference, observed, lcs_graph_gold))
    paths_gold = traversal_gold(reference, observed, lcs_graph_gold, atomics=True)
    hgvs_gold = set()
    for path in paths_gold:
        hgvs_gold.add(to_hgvs(path, reference))
    print(hgvs_gold)
    print(len(hgvs_gold))

    print("test")
    graph = graph_test(reference, observed, nodes_test)
    print(to_dot_test(reference, graph))
    paths_test = traversal(reference, observed, graph, atomics=True)
    hgvs_test = set()
    for path in paths_test:
        hgvs_test.add(to_hgvs(path, reference))
    print(hgvs_test)
    print(len(hgvs_test))

    assert len(list(paths_gold)) == len(list(paths_test))
    assert hgvs_gold == hgvs_test

    print(set(hgvs_test) - set(hgvs_gold))
    print(set(hgvs_gold) - set(hgvs_test))


if __name__ == '__main__':
    main()
