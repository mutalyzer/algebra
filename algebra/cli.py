import argparse
import random
import string
from algebra.utils import hgvs2obs, ref2seq
from algebra.algebra import compare
from algebra.simple import compare as simple_compare


def process(reference, lhs, rhs, debug, check):
    print(f"reference: {reference}")
    print(f"lhs: {lhs}")
    print(f"rhs: {rhs}")

    relation, ops1, ops2 = compare(reference, lhs, rhs, debug)
    print(f"Distance relation: {relation}")

    if check:
        simple_relation, simple_ops1, simple_ops2 = simple_compare(
            reference, lhs, rhs, debug
        )
        print(f"Simple relation: {simple_relation}")

        if ops1 is not None:
            assert simple_ops1 == ops1
        if ops2 is not None:
            assert simple_ops2 == ops2

    return True


def main():

    parser = argparse.ArgumentParser(
        description="A Boolean Algebra for Genetic Variants"
    )

    # Reference
    reference_group = parser.add_mutually_exclusive_group()
    reference_group.add_argument(
        "--reference", type=str, help="reference sequence string", default=None
    )
    reference_group.add_argument(
        "--reference-file", type=str, help="reference sequence from file", default=""
    )
    reference_group.add_argument(
        "--reference-id", type=str, help="reference id", default=""
    )

    # LHS observation
    lhs_group = parser.add_mutually_exclusive_group()
    lhs_group.add_argument(
        "--lhs", type=str, help="lhs observed sequence string", default=None
    )
    lhs_group.add_argument(
        "--lhs-file", type=str, help="lhs observed sequence from file", default=""
    )
    lhs_group.add_argument(
        "--lhs-hgvs", type=str, help="lhs observed sequence hgvs", default=""
    )

    # RHS observation
    second_group = parser.add_mutually_exclusive_group()
    second_group.add_argument(
        "--rhs", type=str, help="rhs observed sequence string", default=None
    )
    second_group.add_argument(
        "--rhs-file", type=str, help="rhs observed sequence from file", default=""
    )
    second_group.add_argument(
        "--rhs-hgvs", type=str, help="rhs observed sequence hgvs", default=""
    )

    # Random sequence options
    parser.add_argument(
        "--alphabet",
        action="store_const",
        const=string.ascii_lowercase,
        default="ACTG",
        required=False,
        help="use the whole alphabet for random strings (default: ACTG)",
    )
    parser.add_argument(
        "--min",
        type=int,
        dest="min_rand",
        help="minimal random string length",
        default="0",
    )
    parser.add_argument(
        "--max",
        type=int,
        dest="max_rand",
        help="maximal random string length",
        default="30",
    )

    # Debugging and validation
    parser.add_argument("--debug", help="print debug output", action="store_true")
    parser.add_argument("--check", help="validate implementations", action="store_true")

    args = parser.parse_args()

    # Use HGVS description
    if args.lhs_hgvs:
        args.lhs = hgvs2obs(args.lhs_hgvs)
    if args.rhs_hgvs:
        args.rhs = hgvs2obs(args.rhs_hgvs)
    if args.reference_id:
        args.reference = ref2seq(args.reference_id)

    # Read sequences from file
    if args.reference_file:
        args.reference = open(args.reference_file).read().strip()
    if args.lhs_file:
        args.lhs = open(args.lhs_file).read().strip()
    if args.rhs_file:
        args.rhs = open(args.rhs_file).read().strip()

    # Use random strings if sequence, description or file is not supplied
    if args.reference is None:
        args.reference = "".join(
            random.choice(args.alphabet)
            for _ in range(random.randint(args.min_rand, args.max_rand))
        )
    if args.lhs is None:
        args.lhs = "".join(
            random.choice(args.alphabet)
            for _ in range(random.randint(args.min_rand, args.max_rand))
        )
    if args.rhs is None:
        args.rhs = "".join(
            random.choice(args.alphabet)
            for _ in range(random.randint(args.min_rand, args.max_rand))
        )

    ret = process(args.reference, args.lhs, args.rhs, args.debug, args.check)

    if not ret:
        exit(-1)
