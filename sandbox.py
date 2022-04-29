import sys
from algebra.utils import random_sequence, random_variants
from algebra.variants.variant import to_hgvs


def main():
    reference = random_sequence(10_000, 1_000)
    variants = list(random_variants(reference))
    print(to_hgvs(variants, reference, sequence_prefix=False))


if __name__ == "__main__":
    main()
