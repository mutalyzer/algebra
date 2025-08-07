from algebra import Variant
from algebra.extractor import extract


def main():
    reference = "GTGTGTTTTTTTAACAGGGA"
    variants = [Variant(4, 5, ""), Variant(12, 13, "")]

    canonical, graph = extract(reference, variants)
    print(canonical)
    print(graph.supremal())
    print(graph.local_supremal())


if __name__ == "__main__":
    main()
