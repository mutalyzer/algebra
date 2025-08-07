from algebra import LCSgraph, Variant


def main():
    graph = LCSgraph.from_variants("AAAAAAA", [Variant(1, 1, "XXX")])
    print(graph)


if __name__ == "__main__":
    main()
