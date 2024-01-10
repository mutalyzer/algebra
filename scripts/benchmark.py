import cProfile
import pstats
import sys
from itertools import combinations
from algebra.lcs import LCSgraph
from algebra.utils import fasta_sequence
from algebra.variants import parse_hgvs, to_hgvs
from algebra.relations.graph_based import compare


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
def graphs(reference, variants):
    for variant in variants:
        graph = LCSgraph.from_variant(reference, variant["variants"])
        variant["graph"] = graph


@benchmark
def pairwise(reference, variants):
    return [{"lhs": lhs["label"],
             "rhs": rhs["label"],
             "relation": compare(reference, lhs["graph"], rhs["graph"])
             } for lhs, rhs in combinations(variants, 2)]


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

    graphs(reference, variants)

    with open("data/benchmark_fast.txt", "w", encoding="utf-8") as file:
        for variant in variants:
            print(variant["label"], f"{ref_seq_id}:g.{to_hgvs(variant['variants'], reference)}", variant["graph"].supremal.to_spdi(reference_id=ref_seq_id), file=file)

    relations = pairwise(reference, variants)
    with open("data/benchmark_fast_relations.txt", "w", encoding="utf-8") as file:
        for entry in relations:
            print(entry["lhs"], entry["rhs"], entry["relation"].value, file=file)


if __name__ == "__main__":
    main()
