from itertools import product
import sys
from algebra.lcs.wupp import edit, lcs_graph, traversal, to_dot
from algebra.relations import Relation, disjoint_variants, ops_set
from algebra.utils import random_sequence


TESTS = [
    ("A", "B", "B", Relation.EQUIVALENT),
    ("AAA", "AAB", "AAB", Relation.EQUIVALENT),
    ("AAA", "AA", "AA", Relation.EQUIVALENT),
    ("AAA", "CAAA", "CAAA", Relation.EQUIVALENT),
    ("AAA", "", "", Relation.EQUIVALENT),
    ("AAA", "ABB", "ABB", Relation.EQUIVALENT),
    ("AA", "AB", "BB", Relation.IS_CONTAINED),
    ("AAA", "AAB", "ABB", Relation.IS_CONTAINED),
    ("", "A", "AA", Relation.IS_CONTAINED),
    ("", "AB", "ABAB", Relation.IS_CONTAINED),
    ("ATATA", "ATATAB", "ATBTAB", Relation.IS_CONTAINED),
    ("", "BB", "BAB", Relation.IS_CONTAINED),
    ("CATATATC", "CATATTATC", "CATATATATC", Relation.IS_CONTAINED),
    ("AA", "BB", "AB", Relation.CONTAINS),
    ("AAA", "ABB", "AAB", Relation.CONTAINS),
    ("", "AA", "A", Relation.CONTAINS),
    ("", "ABAB", "AB", Relation.CONTAINS),
    ("ATATA", "ATBTAB", "ATATAB", Relation.CONTAINS),
    ("", "BAB", "BB", Relation.CONTAINS),
    ("CATATATC", "CATATATATC", "CATATTATC", Relation.CONTAINS),
    ("A", "B", "C", Relation.OVERLAP),
    ("AAA", "ABC", "ABD", Relation.OVERLAP),
    ("AAA", "BBA", "ABB", Relation.OVERLAP),
    ("", "BC", "CAB", Relation.OVERLAP),
    ("ATA", "BTA", "ATB", Relation.DISJOINT),
    ("AAA", "BAA", "AAA", Relation.DISJOINT),
    ("AAA", "AAA", "AAB", Relation.DISJOINT),
    ("AAA", "BAAA", "AAAB", Relation.DISJOINT),
    ("AAA", "AAAB", "BAAA", Relation.DISJOINT),
    ("", "A", "B", Relation.DISJOINT),
    ("T", "GG", "GGTA", Relation.OVERLAP),
    ("TC", "GTC", "GAA", Relation.IS_CONTAINED),
    ("T", "GC", "CT", Relation.CONTAINS),
    ("CT", "TG", "GC", Relation.DISJOINT),
    ("A", "ABD", "ABC", Relation.OVERLAP),
    ("A", "AB", "AC", Relation.DISJOINT),
    ("AAA", "BAAA", "AAAB", Relation.DISJOINT),
    ("A", "BAC", "BAD", Relation.OVERLAP),
    ("AA", "BAAC", "BAAD", Relation.OVERLAP),
    ("AAA", "BAAAC", "BAAAD", Relation.OVERLAP),
    ("TGTA", "CTGCT", "TAGGAACG", Relation.DISJOINT),
    ("CT", "GT", "AT", Relation.OVERLAP),
]


def are_disjoint(lhs, rhs):
    for lhs_variant, rhs_variant in product(lhs, rhs):
        if not disjoint_variants(lhs_variant, rhs_variant):
            return False
    return True


def compare(reference, lhs, rhs):
    print(f'"{reference}" "{lhs}" "{rhs}"')

    lhs_distance, lhs_lcs_nodes = edit(reference, lhs)
    rhs_distance, rhs_lcs_nodes = edit(reference, rhs)

    lhs_root, lhs_edges = lcs_graph(reference, lhs, lhs_lcs_nodes)
    rhs_root, rhs_edges = lcs_graph(reference, rhs, rhs_lcs_nodes)

    print(to_dot(reference, lhs_root))
    print(to_dot(reference, rhs_root))

    print(sorted([variant.to_hgvs(reference) for variant in lhs_edges]))
    print(sorted([variant.to_hgvs(reference) for variant in rhs_edges]))

    lhs_ops = ops_set(lhs_edges)
    rhs_ops = ops_set(rhs_edges)

    print(sorted({variant.to_hgvs(reference) for variant in lhs_ops}))
    print(sorted({variant.to_hgvs(reference) for variant in rhs_ops}))

    print(lhs_ops.isdisjoint(rhs_ops))
    print(are_disjoint(lhs_edges, rhs_edges))

    assert lhs_ops.isdisjoint(rhs_ops) == are_disjoint(lhs_edges, rhs_edges)


def main():
    if len(sys.argv) < 4:
        reference = random_sequence(15)
        lhs = random_sequence(15)
        rhs = random_sequence(15)
    else:
        reference = sys.argv[1]
        lhs = sys.argv[2]
        rhs = sys.argv[3]

    compare(reference, lhs, rhs)

    for reference, lhs, rhs, relation in TESTS:
        compare(reference, lhs, rhs)


if __name__ == "__main__":
    main()
