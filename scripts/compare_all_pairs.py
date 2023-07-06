import argparse
from itertools import combinations
import sys
from algebra import Relation
from algebra.relations.supremal_based import compare, find_supremal
from algebra.utils import fasta_sequence
from algebra.variants import parse_spdi


def main():
    parser = argparse.ArgumentParser(description="Compare all variant pairs from line-based input")
    parser.add_argument("--refseq-id", required=True, type=str, help="RefSeq identifier")
    parser.add_argument("--fasta", required=True, type=str, help="Path to reference fasta file")
    args = parser.parse_args()

    with open(args.fasta, encoding="utf-8") as file:
        reference = fasta_sequence(file)

    variants = []
    for line in sys.stdin:
        variant, *_ = find_supremal(reference, parse_spdi(line.strip())[0])
        variants.append(variant)

    for lhs, rhs in combinations(variants, 2):
        relation = compare(reference, lhs, rhs)
        if relation != Relation.DISJOINT:
            print(lhs.to_spdi(args.refseq_id), rhs.to_spdi(args.refseq_id), relation.value)


if __name__ == "__main__":
    main()
