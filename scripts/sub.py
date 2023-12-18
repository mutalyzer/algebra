from itertools import combinations
from sys import argv
from algebra import Relation, compare
from algebra.extractor import extract, to_hgvs
from algebra.lcs import bfs_traversal, supremal
from algebra.lcs.all_lcs import LCSnode
from algebra.utils import to_dot
from algebra.variants import Variant, parse_hgvs


def unique_matches(graph):
    for node in {sink for _, sink, _ in bfs_traversal(graph)} | {graph}:
        if node.length == 0:
            if not node.edges:
                yield LCSnode(node.row + 1, node.col + 1, 0)
            else:
                yield LCSnode(node.row, node.col, 0)
        else:
            yield from (LCSnode(node.row + i, node.col + i, 1) for i in range(node.length))


def delins(observed, shift, lhs, rhs):
    return Variant(lhs.row + lhs.length, rhs.row + rhs.length - 1,
                   observed[lhs.col + lhs.length - shift:rhs.col + rhs.length - 1 - shift])


def main():
    reference = argv[1]
    minuend = parse_hgvs(argv[2], reference)
    supremal_variant, graph = supremal(reference, minuend)
    print("\n".join(to_dot(reference, graph)))

    seen = set()
    shift = graph.row
    matches = sorted(unique_matches(graph))
    source = matches[0]
    sink = matches[-1]
    for lhs, rhs in combinations(matches, 2):
        if rhs.row > lhs.row and rhs.col > lhs.col:
            variant = delins(supremal_variant.sequence, shift, lhs, rhs)

            subtrahend, *_ = extract(reference, [variant])
            if tuple(subtrahend) in seen:
                continue

            seen.add(tuple(subtrahend))

            difference = []
            if lhs != source:
                prefix = delins(supremal_variant.sequence, shift, source, lhs)
                if prefix:
                    difference.append(prefix)
            if rhs != sink:
                suffix = delins(supremal_variant.sequence, shift, rhs, sink)
                if suffix:
                    difference.append(suffix)

            relation0 = compare(reference, minuend, subtrahend)
            relation1 = compare(reference, minuend, difference)
            relation2 = compare(reference, subtrahend, difference)

            assert relation0 != Relation.IS_CONTAINED
            assert (relation0, relation1) in [
                (Relation.EQUIVALENT, Relation.DISJOINT),
                (Relation.DISJOINT, Relation.EQUIVALENT),
                (Relation.CONTAINS, Relation.CONTAINS),
                (Relation.OVERLAP, Relation.OVERLAP),
            ]
            assert compare(reference, minuend, [variant, *difference]) == Relation.EQUIVALENT
            if relation0 == Relation.CONTAINS:
                assert relation2 == Relation.DISJOINT
            if relation2 == Relation.OVERLAP:
                assert relation0 == Relation.OVERLAP

            print(lhs, rhs, variant.to_hgvs(reference), to_hgvs(subtrahend, reference), relation0, to_hgvs(difference, reference), relation1, relation2)


if __name__ == "__main__":
    main()
