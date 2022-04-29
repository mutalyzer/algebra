import sys
from algebra.lcs.wupp import edit, lcs_graph, to_dot, traversal
from algebra.utils import random_sequence, random_variants
from algebra.variants.variant import patch, to_hgvs


def main():
    reference = random_sequence(10_000, 1_000)
    variants = list(random_variants(reference, p=0.001, mu_deletion=1.5, mu_insertion=1.7))
    observed = patch(reference, variants)

    distance, lcs_nodes = edit(reference, observed)

    print(to_hgvs(variants, reference))
    print(distance)
    print(sum([len(level) for level in lcs_nodes]))

    root, edges = lcs_graph(reference, observed, lcs_nodes)

    print(to_dot(reference, root))

    for description in traversal(root):
        print(to_hgvs(description, reference, sequence_prefix=False))


if __name__ == "__main__":
    main()
