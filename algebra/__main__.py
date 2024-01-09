"""Entry point for the Algebra package.

Provides a command-line interface to interact with sequences and variants
to determine their relations.
"""


import argparse
from algebra.extractor import extract, extract_sequence, local_supremal, to_hgvs as to_hgvs_extractor
from algebra.relations.sequence_based import compare as compare_sequence
from algebra.relations.variant_based import compare
from algebra.utils import fasta_sequence, random_sequence, random_variants, slice_sequence, to_dot
from algebra.variants import parse_hgvs, parse_spdi, patch, to_hgvs


def cli_compare(reference, args):
    """Compare two variants."""
    lhs_is_variant = any([args.lhs_hgvs, args.lhs_spdi, args.lhs_random_variant])
    if args.lhs:
        lhs = args.lhs
    elif args.lhs_hgvs:
        lhs = parse_hgvs(args.lhs_hgvs, reference=reference)
    elif args.lhs_spdi:
        lhs = parse_spdi(args.lhs_spdi)
    elif args.lhs_file:
        with open(args.lhs_file, encoding="utf-8") as file:
            lhs = fasta_sequence(file)
    elif args.lhs_random_variant:
        lhs = list(random_variants(reference, args.random_variant_p))
        print(to_hgvs(lhs, reference))
    else:  # args.lhs_random_sequence:
        lhs = random_sequence(args.random_sequence_max, args.random_sequence_min)
        print(lhs)

    rhs_is_variant = any([args.rhs_hgvs, args.rhs_spdi, args.rhs_random_variant])
    if args.rhs:
        rhs = args.rhs
    elif args.rhs_hgvs:
        rhs = parse_hgvs(args.rhs_hgvs, reference=reference)
    elif args.rhs_spdi:
        rhs = parse_spdi(args.rhs_spdi)
    elif args.rhs_file:
        with open(args.rhs_file, encoding="utf-8") as file:
            rhs = fasta_sequence(file)
    elif args.rhs_random_variant:
        rhs = list(random_variants(reference, args.random_variant_p))
        print(to_hgvs(rhs, reference))
    else:  # args.rhs_random_sequence
        rhs = random_sequence(args.random_sequence_max, args.random_sequence_min)
        print(rhs)

    if lhs_is_variant and rhs_is_variant:
        print(compare(reference, lhs, rhs))
        return

    if lhs_is_variant:
        lhs = patch(reference, lhs)
    elif rhs_is_variant:
        rhs = patch(reference, rhs)

    print(compare_sequence(reference, lhs, rhs))


def cli_extract(reference, args):
    """Extract a canonical variant."""
    is_variant = any([args.observed_hgvs, args.observed_spdi, args.observed_random_variant])
    if args.observed is not None:
        observed = args.observed
    elif args.observed_hgvs:
        observed = parse_hgvs(args.observed_hgvs, reference=reference)
    elif args.observed_spdi:
        observed = parse_spdi(args.observed_spdi)
    elif args.observed_file:
        with open(args.observed_file, encoding="utf-8") as file:
            observed = fasta_sequence(file)
    elif args.observed_random_variant:
        observed = list(random_variants(reference, args.random_variant_p))
        print(to_hgvs(observed, reference))
    else:  # args.observed_random_sequence
        observed = random_sequence(args.random_sequence_max, args.random_sequence_min)
        print(observed)

    if is_variant:
        canonical, graph = extract(reference, observed)
    else:  # observed sequence
        canonical, graph = extract_sequence(reference, observed)

    print(to_hgvs_extractor(canonical, reference))

    if args.all or args.atomics:
        for variants in graph.paths(atomics=args.atomics):
            print(to_hgvs(variants, reference))
    if args.distance:
        print(graph.distance)
    if args.dot:
        print("\n".join(to_dot(reference, graph, atomics=args.atomics)))
    if args.local_supremal:
        print(to_hgvs(local_supremal(reference, graph), reference))
    if args.supremal:
        print(graph.supremal.to_hgvs(reference), graph.supremal.to_spdi(), graph.supremal)


def cli_patch(reference, args):
    """Patch a reference sequence with a variant."""
    if args.hgvs:
        variants = parse_hgvs(args.hgvs, reference=reference)
    elif args.spdi:
        variants = parse_spdi(args.spdi)
    else:  # args.random_variant
        variants = list(random_variants(reference, args.random_variant_p))
        print(to_hgvs(variants, reference))

    print(patch(reference, variants))


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

    extract_parser = commands.add_parser("extract", help="extract a canonical variant")
    extract_parser.add_argument("--all", action="store_true", help="list all minimal variants")
    extract_parser.add_argument("--atomics", action="store_true", help="only deletions and insertions")
    extract_parser.add_argument("--distance", action="store_true", help="output simple edit distance")
    extract_parser.add_argument("--dot", action="store_true", help="output Graphviz DOT")
    extract_parser.add_argument("--local-supremal", action="store_true", help="output local supremal variant")
    extract_parser.add_argument("--supremal", action="store_true", help="output supremal variant")

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

    slice_parser = commands.add_parser("slice", help="slices a reference sequence")

    slice_parser.add_argument("--positions", type=int, nargs="+", required=True, help="positions to slice")
    slice_parser.add_argument("--reverse-complement", action="store_true", help="the reverse complement of the slices")

    args = parser.parse_args()

    if not args.random_sequence_min:
        args.random_sequence_min = args.random_sequence_max

    if args.reference is not None:
        reference = args.reference
    elif args.reference_file:
        with open(args.reference_file, encoding="utf-8") as file:
            reference = fasta_sequence(file)
    else:  # args.reference_random_sequence
        reference = random_sequence(args.random_sequence_max, args.random_sequence_min)
        print(reference)

    if args.command == "compare":
        cli_compare(reference, args)
    elif args.command == "extract":
        cli_extract(reference, args)
    elif args.command == "patch":
        cli_patch(reference, args)
    elif args.command == "slice":
        print(slice_sequence(reference, args.positions, args.reverse_complement))


if __name__ == "__main__":
    main()
