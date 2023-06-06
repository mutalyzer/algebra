from argparse import ArgumentParser

from algebra.lcs import edit, lcs_graph, traversal
from algebra.lcs.all_lcs_dfa import lcs_graph_dfa
from algebra.lcs.all_lcs_dfa import traversal as traversal_no_variant
from algebra.utils import random_sequence


def get_minimal_dfa(reference, observed):
    _, lcs_nodes_dfa = edit(reference, observed)
    source_mdfa, _ = lcs_graph_dfa(reference, observed, lcs_nodes_dfa)

    return traversal_no_variant(source_mdfa)


def get_minimal(reference, observed):
    _, lcs_nodes = edit(reference, observed)
    source, _ = lcs_graph(reference, observed, lcs_nodes)

    return traversal(source)


def compare_minimal(reference, observed):
    minimal_dfa = get_minimal_dfa(reference, observed)
    minimal = get_minimal(reference, observed)

    if sorted(
        minimal_dfa, key=lambda x: (len(x), *((e.start, e.end, e.sequence) for e in x))
    ) != sorted(
        minimal, key=lambda x: (len(x), *((e.start, e.end, e.sequence) for e in x))
    ):
        print("\nminimal dfa:\n")
        print(minimal_dfa)
        print("\n\nminimal regular:\n")
        print(minimal)
        return False
    else:
        return True


def main():
    parser = ArgumentParser(
        "Check and compare the minimal descriptions with the DFA approach"
    )
    parser.add_argument(
        "-s",
        "--sequences",
        action="store",
        nargs=2,
        metavar=("ref", "obs"),
        help="check for the specified sequences",
    )
    parser.add_argument(
        "-i",
        "--interval",
        action="store",
        nargs=2,
        metavar=("max", "min"),
        type=int,
        help="sequences length",
    )

    args = parser.parse_args()
    if args.sequences:
        if compare_minimal(*args.sequences):
            print("equal")
    else:
        count = 0
        len_max = 10 if not args.interval else max(args.interval)
        len_min = 0 if not args.interval else min(args.interval)
        while True:
            count += 1
            ref = random_sequence(len_max, len_min)
            obs = random_sequence(len_max, len_min)
            print(f"\n=== [{len_max}, {len_min}]\n{ref} {obs}\n{len(ref)} {len(obs)}\n---\n{count}")
            assert compare_minimal(ref, obs)


if __name__ == "__main__":
    main()
