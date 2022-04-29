import sys
from algebra.lcs.efficient import edit as edit_gold, build, traversal as traversal_gold
from algebra.lcs.wupp import edit as edit_test, lcs_graph, traversal
from algebra.utils import random_sequence
from algebra.variants.variant import to_hgvs


def main():
    if len(sys.argv) == 1:
        reference = random_sequence(21)
        observed = random_sequence(21)
    elif len(sys.argv) == 3:
        reference = sys.argv[1]
        observed = sys.argv[2]
    else:
        raise("usage")

    print(f'"{reference}" "{observed}"')

    dist_gold, nodes_gold = edit_gold(reference, observed)
    dist_test, nodes_test = edit_test(reference, observed)

    assert dist_test == dist_gold

    _, lcs_graph_gold = build(nodes_gold, reference, observed)
    paths_gold = list(traversal_gold(reference, observed, lcs_graph_gold))
    hgvs_gold = set()
    for variants in paths_gold:
        hgvs_gold.add(to_hgvs(variants, reference))

    atomic_paths_gold = list(traversal_gold(reference, observed, lcs_graph_gold, atomics=True))
    atomic_hgvs_gold = set()
    for variants in atomic_paths_gold:
        atomic_hgvs_gold.add(to_hgvs(variants, reference))

    root, _ = lcs_graph(reference, observed, nodes_test)
    paths_test = list(traversal(root))
    hgvs_test = set()
    for variants in paths_test:
        hgvs_test.add(to_hgvs(variants, reference))

    atomic_paths_test = list(traversal(root, atomics=True))
    atomic_hgvs_test = set()
    for variants in atomic_paths_test:
        atomic_hgvs_test.add(to_hgvs(variants, reference))

    assert len(paths_gold) == len(paths_test)
    assert len(atomic_paths_gold) == len(atomic_paths_test)
    assert hgvs_gold == hgvs_test
    assert atomic_hgvs_gold == atomic_hgvs_test


if __name__ == '__main__':
    main()
