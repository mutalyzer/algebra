import cProfile
import pstats
import sys
from itertools import combinations, product
from algebra import Relation
from algebra.lcs import bfs_traversal, dfs_traversal, edit, supremal
from algebra.utils import fasta_sequence
from algebra.variants import parse_hgvs, to_hgvs


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
def supremals(reference, variants):
    for variant in variants:
        supremal_variant, graph = supremal(reference, variant["variants"])
        variant["supremal"] = supremal_variant
        variant["graph"] = graph


@benchmark
def pairwise(reference, variants):
    return [{"lhs": lhs["label"], "rhs": rhs["label"], "relation": compare_graph(reference, lhs["supremal"], lhs["graph"], rhs["supremal"], rhs["graph"])} for lhs, rhs in combinations(variants, 2)]


def main():
    ref_seq_id = "NC_000022.11"
    with open(f"data/{ref_seq_id}.fasta", encoding="utf-8") as file:
        reference = fasta_sequence(file)

    variants = []
    with open("data/benchmark.txt", encoding="utf-8") as file:
        for line in file:
            label, hgvs = line.split()
            variant = parse_hgvs(hgvs, reference)
            variants.append({
                "label": label,
                "variants": variant,
            })

    supremals(reference, variants)

    with open("data/benchmark_fast.txt", "w", encoding="utf-8") as file:
        for variant in variants:
            print(variant["label"], f"{ref_seq_id}:g.{to_hgvs(variant['variants'], reference)}", variant["supremal"].to_spdi(reference_id=ref_seq_id), file=file)

    relations = pairwise(reference, variants)
    with open("data/benchmark_fast_relations.txt", "w", encoding="utf-8") as file:
        for entry in relations:
            print(entry["lhs"], entry["rhs"], entry["relation"].value, file=file)


if __name__ == "__main__":
    main()
