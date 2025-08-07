from algebra import LCSgraph, Variant


def main():
    reference = "GTGTGTTTTTTTAACAGGGA"
    variants = [Variant(4, 5, ""), Variant(12, 13, "")]

    graph = LCSgraph.from_variants(reference, variants)

    print(graph.supremal())
    print(graph.local_supremal())
    print(graph.canonical())


if __name__ == "__main__":
    main()
