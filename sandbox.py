from algebra.lcs.efficient import edit, build


def main():
    reference = "CATATATCG"
    observed = "CTTATAGCAT"

    distance, lcs_nodes = edit(reference, observed)
    _, graph = build(lcs_nodes, reference, observed)

    for level in graph:
        print([str(node) for node in level])


if __name__ == "__main__":
    main()

