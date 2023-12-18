from itertools import combinations
from sys import argv
from algebra.lcs import bfs_traversal, supremal
from algebra.utils import to_dot
from algebra.variants import Variant, parse_hgvs


def unique_matches(graph):
    for node in {sink for _, sink, _ in bfs_traversal(graph)} | {graph}:
        if node.length == 0:
            if node == graph:
                yield node.row - 1, node.col - 1
            elif not node.edges:
                yield node.row, node.col + 1
            else:
                raise ValueError
        for i in range(node.length):
            yield node.row + i, node.col + i


def main():
    reference = argv[1]
    variants = parse_hgvs(argv[2], reference)
    supremal_variant, graph = supremal(reference, variants)
    print("\n".join(to_dot(reference, graph)))

    shift = graph.row
    for lhs, rhs in combinations(sorted(unique_matches(graph)), 2):
        if rhs[0] > lhs[0] and rhs[1] > lhs[1]:
            print(lhs, rhs)
            variant = Variant(lhs[0] + 1, rhs[0], supremal_variant.sequence[lhs[1] + 1 - shift:rhs[1] - shift])
            print(variant.to_hgvs(reference))


if __name__ == "__main__":
    main()
