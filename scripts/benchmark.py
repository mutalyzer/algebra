import cProfile
import pstats
import sys
from itertools import combinations, product
from algebra import Relation
from algebra.lcs import bfs_traversal, dfs_traversal, edit, supremal
from algebra.utils import fasta_sequence
from algebra.variants import parse_hgvs, to_hgvs


def compare_graph(reference, lhs, rhs):
    if lhs.supremal == rhs.supremal:
        return Relation.EQUIVALENT

    if lhs.supremal.is_disjoint(rhs.supremal):
        return Relation.DISJOINT

    start = min(lhs.supremal.start, rhs.supremal.start)
    end = max(lhs.supremal.end, rhs.supremal.end)

    # TODO: class method on LCSgraph?
    def observed(supremal):
        return (reference[min(start, supremal.start):supremal.start] +
                supremal.sequence +
                reference[supremal.end:max(end, supremal.end)])

    lhs_observed = observed(lhs.supremal)
    rhs_observed = observed(rhs.supremal)
    distance = edit(lhs_observed, rhs_observed)

    if lhs.distance + rhs.distance == distance:
        return Relation.DISJOINT

    if lhs.distance - rhs.distance == distance:
        return Relation.CONTAINS

    if rhs.distance - lhs.distance == distance:
        return Relation.IS_CONTAINED

    for lhs_variant, rhs_variant in product(lhs.edges(), rhs.edges()):
        if not lhs_variant.is_disjoint(rhs_variant):
            return Relation.OVERLAP

    return Relation.DISJOINT


BENCHMARK_ENABLE = False
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
def pairwise(reference, graphs):
    return [{"lhs": lhs.label, "rhs": rhs.label, "relation": compare_graph(reference, lhs, rhs)} for lhs, rhs in combinations(graphs, 2)]


class LCSgraph:
    def __init__(self, label, variants, supremal, graph):
        self.label = label
        self.variants = variants
        self.supremal = supremal
        self.graph = graph

    def edges(self):
        return {edge[0] for *_, edge in bfs_traversal(self.graph)}

    @property
    def distance(self):
        return sum(len(variant) for variant in next(dfs_traversal(self.graph)))


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

    graphs = [LCSgraph(variant["label"], variant["variants"], variant["supremal"], variant["graph"]) for variant in variants]

    with open("data/benchmark_fast.txt", "w", encoding="utf-8") as file:
        for graph in graphs:
            print(graph.label, f"{ref_seq_id}:g.{to_hgvs(graph.variants, reference)}", graph.supremal.to_spdi(reference_id=ref_seq_id), file=file)

    relations = pairwise(reference, graphs)
    with open("data/benchmark_fast_relations.txt", "w", encoding="utf-8") as file:
        for entry in relations:
            print(entry["lhs"], entry["rhs"], entry["relation"].value, file=file)


if __name__ == "__main__":
    main()
