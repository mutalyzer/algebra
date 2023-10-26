from argparse import ArgumentParser
from sys import stdin
from algebra import Variant
from algebra.extractor import extract_supremal, to_hgvs
from algebra.lcs import supremal
from algebra.utils import fasta_sequence
from algebra.variants import parse_spdi, to_hgvs as to_hgvs_simple


def main():
    parser = ArgumentParser(description="Extract variants from line-based input")
    parser.add_argument("--fasta", required=True, type=str, help="Path to reference fasta file")
    args = parser.parse_args()

    with open(args.fasta, encoding="utf-8") as file:
        reference = fasta_sequence(file)

    for line in stdin:
        spdi = line.strip()
        variant = parse_spdi(spdi)[0]
        print(spdi, variant, end=" ")
        supremal_variant, *_ = supremal(reference, variant)
        print(supremal_variant, end=" ")
        canonical = extract_supremal(reference, supremal_variant)
        print(to_hgvs_simple([variant], reference), to_hgvs(canonical, reference))


if __name__ == "__main__":
    main()
