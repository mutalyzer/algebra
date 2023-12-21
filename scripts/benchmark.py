import cProfile
import pstats
import sys
from itertools import combinations, product
from algebra import Relation
from algebra.lcs import bfs_traversal, dfs_traversal, edit, supremal
from algebra.utils import fasta_sequence
from algebra.variants import parse_hgvs


def compare_graph(reference, lhs, lhs_graph, rhs, rhs_graph):
    if lhs == rhs:
        return Relation.EQUIVALENT

    if lhs.is_disjoint(rhs):
        return Relation.DISJOINT

    lhs_distance = sum(len(variant) for variant in next(dfs_traversal(lhs_graph), []))
    rhs_distance = sum(len(variant) for variant in next(dfs_traversal(rhs_graph), []))

    start = min(lhs.start, rhs.start)
    end = max(lhs.end, rhs.end)
    lhs_observed = reference[min(start, lhs.start):lhs.start] + lhs.sequence + reference[lhs.end:max(end, lhs.end)]
    rhs_observed = reference[min(start, rhs.start):rhs.start] + rhs.sequence + reference[rhs.end:max(end, rhs.end)]
    distance = edit(lhs_observed, rhs_observed)

    if lhs_distance + rhs_distance == distance:
        return Relation.DISJOINT

    if lhs_distance - rhs_distance == distance:
        return Relation.CONTAINS

    if rhs_distance - lhs_distance == distance:
        return Relation.IS_CONTAINED

    lhs_edges = {edge[0] for *_, edge in bfs_traversal(lhs_graph)}
    rhs_edges = {edge[0] for *_, edge in bfs_traversal(rhs_graph)}

    for lhs_variant, rhs_variant in product(lhs_edges, rhs_edges):
        if not lhs_variant.is_disjoint(rhs_variant):
            return Relation.OVERLAP

    return Relation.DISJOINT


BENCHMARK_ENABLE = True
BENCHMARK_STATS = "tottime"


def benchmark(func):
    def wrap(*args, **kwargs):
        if BENCHMARK_ENABLE:
            print(f"PROFILE {func.__name__}", file=sys.stderr)
            pr = cProfile.Profile()
            pr.enable()
            result = func(*args, **kwargs)
            pr.disable()
            ps = pstats.Stats(pr, stream=sys.stderr).sort_stats(BENCHMARK_STATS)
            ps.print_stats()
            return result
        return func(*args, **kwargs)
    return wrap


@benchmark
def parsing(reference, raw):
    return [parse_hgvs(variants, reference) for variants in raw.values()]


@benchmark
def supremals(reference, variants):
    return [supremal(reference, variant) for variant in variants]


@benchmark
def pairwise(reference, supremal_variants):
    return [compare_graph(reference, lhs_supremal, lhs_graph, rhs_supremal, rhs_graph) for (lhs_supremal, lhs_graph), (rhs_supremal, rhs_graph) in combinations(supremal_variants, 2)]


def main():
    with open("data/NC_000022.11.fasta", encoding="utf-8") as file:
        reference = fasta_sequence(file)

    with open("data/benchmark.txt", encoding="utf-8") as file:
        raw = {line.split()[0]: line.split()[1] for line in file}

    variants = parsing(reference, raw)
    with open("data/benchmark_variants.txt", "w", encoding="utf-8") as file:
        for variant in variants:
            print(variant, file=file)

    supremal_variants = supremals(reference, variants)
    with open("data/benchmark_supremals.txt", "w", encoding="utf-8") as file:
        for variant in supremal_variants:
            print(variant[0], file=file)

    relations = pairwise(reference, supremal_variants)
    with open("data/benchmark_relations.txt", "w", encoding="utf-8") as file:
        for relation in relations:
            print(relation, file=file)


if __name__ == "__main__":
    main()
