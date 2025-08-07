from algebra import LCSgraph, Variant


def main():
    reference = "GTGTGTTTTTTTAACAGGGA"
    variants = [Variant(8, 9, "")]

    graph = LCSgraph.from_variants(reference, variants)

    print(graph.supremal())
    print(graph.local_supremal())


if __name__ == "__main__":
    main()
