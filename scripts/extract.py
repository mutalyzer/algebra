from argparse import ArgumentParser
from sys import stdin
from algebra.extractor import extract_supremal
from algebra.relations.supremal_based import find_supremal
from algebra.utils import fasta_sequence
from algebra.variants import parse_spdi


def main():
    parser = ArgumentParser(description="Extract variants from line-based input")
    parser.add_argument("--refseq-id", required=True, type=str, help="RefSeq identifier")
    parser.add_argument("--fasta", required=True, type=str, help="Path to reference fasta file")
    args = parser.parse_args()

    with open(args.fasta, encoding="utf-8") as file:
        reference = fasta_sequence(file)

    for line in stdin:
        spdi = line.strip()
        variant = parse_spdi(spdi)[0]
        supremal = find_supremal(reference, variant)
        canonical = list(extract_supremal(reference, supremal))

        print(spdi, variant, supremal, canonical)
        return


if __name__ == "__main__":
    main()
