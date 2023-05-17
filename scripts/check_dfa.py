from argparse import ArgumentParser

from algebra.lcs import edit, lcs_graph, traversal
from algebra.lcs.all_lcs_dfa import edit as edit_mdfa
from algebra.lcs.all_lcs_dfa import lcs_graph_dfa
from algebra.lcs.all_lcs_dfa import traversal as traversal_no_variant
from algebra.utils import random_sequence
from algebra.variants.variant import to_hgvs


def get_minimal_mdfa(reference, observed):
    _, lcs_nodes_mdfa = edit_mdfa(reference, observed)
    source_mdfa, sink_mdfa = lcs_graph_dfa(reference, observed, lcs_nodes_mdfa)
    return traversal_no_variant(source_mdfa)


def get_minimal(reference, observed):
    _, lcs_nodes = edit(reference, observed)
    source, sink = lcs_graph(reference, observed, lcs_nodes)

    return traversal(source)


def compare_minimal(reference, observed):
    minimal_dfa = sorted(
        [to_hgvs(vars, reference) for vars in get_minimal_mdfa(reference, observed)],
        key=str,
    )
    minimal = sorted(
        [to_hgvs(vars, reference) for vars in get_minimal(reference, observed)], key=str
    )
    if minimal != minimal_dfa:
        print("\nminimal dfa:\n")
        print(minimal_dfa)
        print("\n\nminimal regular:\n")
        print(minimal)
        print("\n\n minimal - minimal_dfa:")
        print(set(minimal) - set(minimal_dfa))
        print("\n\n minimal_dfa - minimal:")
        print(set(minimal_dfa) - set(minimal))
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

    args = parser.parse_args()
    if args.sequences:
        if compare_minimal(*args.sequences):
            print("equal")
    else:
        count = 0
        while True:
            count += 1
            ref = random_sequence(25, 10)
            obs = random_sequence(25, 10)
            print(f"\n===\n{ref} {obs}\n---\n{count}")
            assert compare_minimal(ref, obs)


if __name__ == "__main__":
    main()
