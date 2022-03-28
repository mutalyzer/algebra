from algebra.lcs.efficient import edit, lcs_graph


def main():
    reference = "CATATATCG"
    observed = "CTTATAGCAT"

    distance, lcs_nodes = edit(reference, observed)
    print(lcs_nodes)
    graph = lcs_graph(reference, observed, lcs_nodes)


if __name__ == "__main__":
    main()

