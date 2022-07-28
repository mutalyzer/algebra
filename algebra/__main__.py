"""Entry point for the Algebra package.

Provides a command-line interface to interact with variants and determine
their relations.
"""


import argparse
from .lcs import edit, lcs_graph, traversal
from .relations.sequence_based import compare
from .relations.supremal_based import spanning_variant
from .utils import fasta_sequence, random_sequence, random_variants, to_dot
from .variants import parse_hgvs, parse_spdi, patch, to_hgvs


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(description="A Boolean Algebra for Genetic Variants")
    parser.add_argument("--random-sequence-min", type=int, help="minimum length for random sequences")
    parser.add_argument("--random-sequence-max", type=int, default=1_000, help="maximum length for random sequences")

    parser.add_argument("--random-variant-p", type=float, help="change per base of a variant")

    reference_group = parser.add_mutually_exclusive_group()
    reference_group.add_argument("--reference", type=str, help="a reference sequence as string")
    reference_group.add_argument("--reference-file", type=str, help="a reference sequence from a file")
    reference_group.add_argument("--reference-random-sequence", action="store_true", help="a random reference sequence (default)")

    # https://bugs.python.org/issue29298
    commands = parser.add_subparsers(dest="command", required=True, help="Commands")

    compare_parser = commands.add_parser("compare", help="compare two variants")

    lhs_group = compare_parser.add_mutually_exclusive_group()
    lhs_group.add_argument("--lhs", type=str, help="an observed sequence as string (lhs)")
    lhs_group.add_argument("--lhs-hgvs", type=str, help="a variant in HGVS (lhs)")
    lhs_group.add_argument("--lhs-spdi", type=str, help="a variant in SPDI (lhs)")
    lhs_group.add_argument("--lhs-file", type=str, help="an observed sequence from a file (lhs)")
    lhs_group.add_argument("--lhs-random-variant", action="store_true", help="a random variant")
    lhs_group.add_argument("--lhs-random-sequence", action="store_true", help="a random sequence (default)")

    rhs_group = compare_parser.add_mutually_exclusive_group()
    rhs_group.add_argument("--rhs", type=str, help="an observed sequence as string (rhs)")
    rhs_group.add_argument("--rhs-hgvs", type=str, help="a variant in HGVS (rhs)")
    rhs_group.add_argument("--rhs-spdi", type=str, help="a variant in SPDI (rhs)")
    rhs_group.add_argument("--rhs-file", type=str, help="an observed sequence from a file (rhs)")
    rhs_group.add_argument("--rhs-random-variant", action="store_true", help="a random variant")
    rhs_group.add_argument("--rhs-random-sequence", action="store_true", help="a random sequence (default)")

    extract_parser = commands.add_parser("extract", help="extract an HGVS description from an observed sequence")
    extract_parser.add_argument("--all", action="store_true", help="list all minimal HGVS descriptions (default)")
    extract_parser.add_argument("--atomics", action="store_true", help="only deletions and insertions")
    extract_parser.add_argument("--distance", action="store_true", help="output simple edit distance")
    extract_parser.add_argument("--dot", action="store_true", help="output Graphviz DOT")
    extract_parser.add_argument("--supremal-variant", action="store_true", help="output supremal variant")

    observed_group = extract_parser.add_mutually_exclusive_group()
    observed_group.add_argument("--observed", type=str, help="an observed sequence as string")
    observed_group.add_argument("--observed-hgvs", type=str, help="a variant in HGVS")
    observed_group.add_argument("--observed-spdi", type=str, help="a variant in SPDI")
    observed_group.add_argument("--observed-file", type=str, help="an observed sequence from a file")
    observed_group.add_argument("--observed-random-variant", action="store_true", help="a random variant")
    observed_group.add_argument("--observed-random-sequence", action="store_true", help="a random sequence (default)")

    patch_parser = commands.add_parser("patch", help="patch a reference sequence with a variant")

    variant_group = patch_parser.add_mutually_exclusive_group()
    variant_group.add_argument("--hgvs", type=str, help="a variant in HGVS")
    variant_group.add_argument("--spdi", type=str, help="a variant in SPDI")
    variant_group.add_argument("--random-variant", action="store_true", help="a random variant (default)")

    args = parser.parse_args()

    if not args.random_sequence_min:
        args.random_sequence_min = args.random_sequence_max

    if args.reference is not None:
        reference = args.reference
    elif args.reference_file is not None:
        with open(args.reference_file, encoding="utf-8") as file:
            reference = fasta_sequence(file.readlines())
    else:  # args.reference_random_sequence
        reference = random_sequence(args.random_sequence_max, args.random_sequence_min)
        print(reference)

    if args.command == "compare":
        if args.lhs is not None:
            lhs = args.lhs
        elif args.lhs_hgvs is not None:
            lhs = patch(reference, parse_hgvs(args.lhs_hgvs, reference=reference))
        elif args.lhs_spdi is not None:
            lhs = patch(reference, parse_spdi(args.lhs_spdi))
        elif args.lhs_file is not None:
            with open(args.lhs_file, encoding="utf-8") as file:
                lhs = fasta_sequence(file.readlines())
        elif args.lhs_random_variant:
            variants = list(random_variants(reference, args.random_variant_p))
            lhs = patch(reference, variants)
            print(to_hgvs(variants, reference))
        else:  # args.lhs_random_sequence:
            lhs = random_sequence(args.random_sequence_max, args.random_sequence_min)
            print(lhs)

        if args.rhs is not None:
            rhs = args.rhs
        elif args.rhs_hgvs is not None:
            rhs = patch(reference, parse_hgvs(args.rhs_hgvs, reference=reference))
        elif args.rhs_spdi is not None:
            rhs = patch(reference, parse_spdi(args.rhs_spdi))
        elif args.rhs_file is not None:
            with open(args.rhs_file, encoding="utf-8") as file:
                rhs = fasta_sequence(file.readlines())
        elif args.rhs_random_variant:
            variants = list(random_variants(reference, args.random_variant_p))
            rhs = patch(reference, variants)
            print(to_hgvs(variants, reference))
        else:  # args.rhs_random_sequence
            rhs = random_sequence(args.random_sequence_max, args.random_sequence_min)
            print(rhs)

        print(compare(reference, lhs, rhs))

    elif args.command == "extract":
        if args.observed is not None:
            observed = args.observed
        elif args.observed_hgvs is not None:
            observed = patch(reference, parse_hgvs(args.observed_hgvs, reference=reference))
        elif args.observed_spdi is not None:
            observed = patch(reference, parse_spdi(args.observed_spdi))
        elif args.observed_file is not None:
            with open(args.observed_file, encoding="utf-8") as file:
                observed = fasta_sequence(file.readlines())
        elif args.observed_random_variant:
            variants = list(random_variants(reference, args.random_variant_p))
            observed = patch(reference, variants)
            print(to_hgvs(variants, reference))
        else:  # args.observed_random_sequence
            observed = random_sequence(args.random_sequence_max, args.random_sequence_min)
            print(observed)

        distance, lcs_nodes = edit(reference, observed)
        root, edges = lcs_graph(reference, observed, lcs_nodes)

        if args.distance:
            print(distance)
        if args.dot:
            print(to_dot(reference, root))
        if args.supremal_variant:
            variant = spanning_variant(reference, observed, edges)
            print(variant.to_hgvs(reference), variant)
        if True or args.all:
            for variants in traversal(root, args.atomics):
                print(to_hgvs(variants, reference))

    elif args.command == "patch":
        if args.hgvs is not None:
            variants = parse_hgvs(args.hgvs, reference=reference)
        elif args.spdi is not None:
            variants = parse_spdi(args.spdi)
        else:  # args.random_variant
            variants = list(random_variants(reference, args.random_variant_p))
            print(to_hgvs(variants, reference))

        print(patch(reference, variants))


if __name__ == "__main__":
    main()
