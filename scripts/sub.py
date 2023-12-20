from itertools import combinations
from sys import argv
from algebra import Relation, compare
from algebra.extractor import extract, to_hgvs
from algebra.lcs import bfs_traversal, supremal
from algebra.lcs.all_lcs import LCSnode
from algebra.utils import to_dot
from algebra.variants import Variant, parse_hgvs, to_hgvs as to_hgvs_simple


def unique_matches(graph):
    for node in {sink for _, sink, _ in bfs_traversal(graph)} | {graph}:
        if node == graph:
            yield LCSnode(node.row, node.col, 0)
        if not node.edges:
            yield LCSnode(node.row + node.length, node.col + node.length, 0)
        yield from (LCSnode(node.row + i, node.col + i, 1) for i in range(node.length))


def delins(observed, shift, lhs, rhs):
    return Variant(lhs.row + lhs.length, rhs.row, observed[lhs.col + lhs.length - shift:rhs.col - shift])


def subtract(reference, minuend):

    supremal_variant, graph = supremal(reference, minuend)
    print("\n".join(to_dot(reference, graph, labels=False)))

    seen = set()
    shift = graph.row
    matches = sorted(unique_matches(graph))
    source = matches[0]
    sink = matches[-1]
    for lhs, rhs in combinations(matches, 2):
        if rhs.row < lhs.row + lhs.length or rhs.col < lhs.col + lhs.length:
            continue

        variant = delins(supremal_variant.sequence, shift, lhs, rhs)

        subtrahend, *_ = extract(reference, [variant])
        if tuple(subtrahend) in seen:
            continue
        seen.add(tuple(subtrahend))

        difference = []
        if lhs.length:
            difference.append(delins(supremal_variant.sequence, shift, source, lhs))
        if rhs.length:
            difference.append(delins(supremal_variant.sequence, shift, rhs, sink))
        difference_norm, *_ = extract(reference, difference)

        assert compare(reference, minuend, [variant, *difference]) == Relation.EQUIVALENT

        relation0 = compare(reference, minuend, subtrahend)
        relation1 = compare(reference, minuend, difference)
        relation2 = compare(reference, subtrahend, difference)

        print(lhs, rhs, variant, variant.to_hgvs(reference), to_hgvs(subtrahend, reference), relation0, to_hgvs(difference_norm, reference), to_hgvs_simple(difference, reference), relation1, relation2)


def main():
    reference = argv[1]
    minuend = parse_hgvs(argv[2], reference)

    subtract(reference, minuend)


if __name__ == "__main__":
    main()
