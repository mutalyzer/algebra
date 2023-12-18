from itertools import combinations
from sys import argv
from algebra.lcs import bfs_traversal, supremal
from algebra.variants import Variant, parse_hgvs


def unique_matches(graph):
    results = []
    for node in {sink for _, sink, _ in bfs_traversal(graph)} | {graph}:
        for i in range(node.length):
            results.append((node.row + i, node.col + i))
        if node.length == 0:
            assert node == graph or not node.edges
            results.append((node.row, node.col))

    return sorted(results)


def dominates(lhs, rhs):
    return lhs[0] > rhs[0] and lhs[1] > rhs[1]


def delins(observed, lhs, rhs):
    return Variant(lhs[0] + 1, rhs[0], observed[lhs[1]:rhs[1] - 1])


def main():
    reference = argv[1]
    variants = parse_hgvs(argv[2], reference)

    sup, graph = supremal(reference, variants)

    for lhs, rhs in combinations(unique_matches(graph), 2):

        if dominates(lhs, rhs):
            lhs, rhs = rhs, lhs

        if dominates(rhs, lhs):
            variant = delins(sup.sequence, lhs, rhs)
            print(lhs, rhs)
            print(variant.to_hgvs(reference))


if __name__ == "__main__":
    main()
