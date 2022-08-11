from argparse import ArgumentParser
from sys import stdin
from algebra import Variant
from algebra.extractor import extract_supremal, to_hgvs
from algebra.relations.supremal_based import find_supremal
from algebra.utils import fasta_sequence
from algebra.variants import parse_spdi, to_hgvs as to_hgvs_simple


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
        print(spdi, variant, end=" ", flush=True)
        to_skip = [
            Variant(117642410, 117643194, "AGT"), # NC
            Variant(181627, 182411, "AGT"), # NG
            Variant(117664686, 117665566, "GGT"), # NC
            Variant(203903, 204783, "GGT"), # NG
        ]
        if variant in to_skip:
            print("skipped")
            continue
        supremal = find_supremal(reference, variant)
        print(supremal, end=" ", flush=True)
        canonical = list(extract_supremal(reference, supremal))
        print(to_hgvs_simple([variant], reference), to_hgvs(canonical, reference))


if __name__ == "__main__":
    main()
