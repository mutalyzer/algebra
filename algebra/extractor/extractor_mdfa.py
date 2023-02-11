from argparse import ArgumentParser

from algebra.lcs.all_lcs_mdfa import edit, lcs_graph_mdfa, to_dot


def extract(reference, observed):
    _, lcs_nodes = edit(reference, observed)
    source, sink = lcs_graph_mdfa(reference, observed, lcs_nodes)

    return to_dot(reference, source, sink)


def main():
    # print(extract("TTT", "TTTTAT"))
    parser = ArgumentParser(description="Extract variants from line-based input")
    parser.add_argument("ref", type=str, help="reference sequence")
    parser.add_argument("obs", type=str, help="observed sequence")

    args = parser.parse_args()

    print(extract(args.ref, args.obs))


if __name__ == "__main__":
    main()
