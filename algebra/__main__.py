import argparse
from algebra.lcs.all_lcs import edit, lcs_graph, traversal, to_dot
from algebra.relations import compare
from algebra.utils import random_sequence, random_variants
from algebra.variants.parser import Parser
from algebra.variants.variant import patch, to_hgvs


def main():
    parser = argparse.ArgumentParser(description="A Boolean Algebra for Genetic Variants")
    parser.add_argument("--random-min", type=int, help="minimum length for random sequences")
    parser.add_argument("--random-max", type=int, default=1_000, help="maximum length for random sequences")
    parser.add_argument("--random-p", type=float, help="change per base of a variant")

    reference_group = parser.add_mutually_exclusive_group(required=True)
    reference_group.add_argument("--reference", type=str, help="a reference sequence as string")
    reference_group.add_argument("--reference-file", type=str, help="a reference sequence from a file")
    reference_group.add_argument("--reference-random", action="store_true", help="a random reference sequence")

    # https://bugs.python.org/issue29298
    commands = parser.add_subparsers(dest="command", required=True, help="Commands")

    compare_parser = commands.add_parser("compare", help="compare two variants")

    lhs_group = compare_parser.add_mutually_exclusive_group(required=True)
    lhs_group.add_argument("--lhs", type=str, help="an observed sequence as string (lhs)")
    lhs_group.add_argument("--lhs-hgvs", type=str, help="an HGVS variant (lhs)")
    lhs_group.add_argument("--lhs-spdi", type=str, help="a SPDI variant (lhs)")
    lhs_group.add_argument("--lhs-file", type=str, help="an observed sequence from a file (lhs)")
    lhs_group.add_argument("--lhs-random", action="store_true", help="a random variant")

    rhs_group = compare_parser.add_mutually_exclusive_group(required=True)
    rhs_group.add_argument("--rhs", type=str, help="an observed sequence as string (rhs)")
    rhs_group.add_argument("--rhs-hgvs", type=str, help="an HGVS variant (rhs)")
    rhs_group.add_argument("--rhs-spdi", type=str, help="a SPDI variant (rhs)")
    rhs_group.add_argument("--rhs-file", type=str, help="an observed sequence from a file (rhs)")
    rhs_group.add_argument("--rhs-random", action="store_true", help="a random variant")

    extract_parser = commands.add_parser("extract", help="extract an HGVS description from an observed sequence")
    extract_parser.add_argument("--all", action="store_true", help="list all minimal HGVS descriptions")
    extract_parser.add_argument("--atomics", action="store_true", help="only deletions and insertions")
    extract_parser.add_argument("--distance", action="store_true", help="output simple edit distance")
    extract_parser.add_argument("--dot", action="store_true", help="output Graphviz DOT")

    observed_group = extract_parser.add_mutually_exclusive_group(required=True)
    observed_group.add_argument("--observed", type=str, help="an observed sequence as string")
    observed_group.add_argument("--observed-hgvs", type=str, help="an HGVS variant")
    observed_group.add_argument("--observed-spdi", type=str, help="a SPDI variant")
    observed_group.add_argument("--observed-file", type=str, help="an observed sequence from a file")
    observed_group.add_argument("--observed-random", action="store_true", help="a random observed sequence")

    patch_parser = commands.add_parser("patch", help="patch a reference sequence with a variant")

    variant_group = patch_parser.add_mutually_exclusive_group(required=True)
    variant_group.add_argument("--hgvs", type=str, help="a variant in HGVS")
    variant_group.add_argument("--spdi", type=str, help="a variant in SPDI")
    variant_group.add_argument("--random-variant", action="store_true", help="a random variant")

    args = parser.parse_args()

    if not args.random_min:
        args.random_min = args.random_max

    if args.reference:
        reference = args.reference
    elif args.reference_file:
        with open(args.reference_file, encoding="utf-8") as file:
            reference = file.read().strip()
    elif args.reference_random:
        reference = random_sequence(args.random_max, args.random_min)
        print(reference)

    if args.command == "compare":
        if args.lhs:
            lhs = args.lhs
        elif args.lhs_hgvs:
            lhs = Parser(args.lhs_hgvs).hgvs()
        elif args.lhs_spdi:
            lhs = Parser(args.lhs_spdi).spdi()
        elif args.lhs_file:
            with open(args.lhs_file, encoding="utf-8") as file:
                lhs = file.read().strip()
        elif args.lhs_random:
            variants = list(random_variants(reference, args.random_p))
            print(to_hgvs(variants, reference, sequence_prefix=False))
            lhs = patch(reference, variants)

        if args.rhs:
            rhs = args.rhs
        elif args.rhs_hgvs:
            rhs = Parser(args.rhs_hgvs).hgvs()
        elif args.rhs_spdi:
            rhs = Parser(args.rhs_spdi).spdi()
        elif args.rhs_file:
            with open(args.rhs_file, encoding="utf-8") as file:
                rhs = file.read().strip()
        elif args.rhs_random:
            variants = list(random_variants(reference, args.random_p))
            print(to_hgvs(variants, reference, sequence_prefix=False))
            rhs = patch(reference, variants)

        print(compare(reference, lhs, rhs))

    elif args.command == "extract":
        if args.observed:
            observed = args.observed
        elif args.observed_hgvs:
            variants = Parser(args.observed_hgvs).hgvs()
            observed = patch(reference, variants)
        elif args.observed_spdi:
            variants = Parser(args.observed_spdi).spdi()
            observed = patch(reference, variants)
        elif args.observed_file:
            with open(args.observed_file, encoding="utf-8") as file:
                observed = file.read().strip()
        elif args.observed_random:
            observed = random_sequence(args.random_max, args.random_min)
            print(observed)

        distance, lcs_nodes = edit(reference, observed)
        root, _ = lcs_graph(reference, observed, lcs_nodes)

        if args.distance:
            print(distance)
        if args.dot:
            print(to_dot(reference, root))
        if True or args.all:
            for variants in traversal(root, args.atomics):
                print(to_hgvs(variants, reference, sequence_prefix=False))

    elif args.command == "patch":
        if args.hgvs:
            variants = Parser(args.hgvs).hgvs()
        elif args.spdi:
            variants = Parser(args.spdi).spdi()
        elif args.random_variant:
            variants = list(random_variants(reference, args.random_p))
            print(to_hgvs(variants, reference, sequence_prefix=False))

        print(patch(reference, variants))


if __name__ == "__main__":
    main()
