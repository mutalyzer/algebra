from algebra.lcs.efficient import edit, build, traversal
from algebra.variants.variant import Variant, to_hgvs
import sys


def main():
    if len(sys.argv) < 3:
        reference = "CATATATCG"
        observed = "CTTATAGCAT"
    else:
        reference = sys.argv[1]
        observed = sys.argv[2]

    distance, lcs_nodes = edit(reference, observed)
    _, graph = build(lcs_nodes, reference, observed)

    for path in sorted(traversal(reference, observed, graph, atomics=True)):
        print(to_hgvs(path, reference))


if __name__ == "__main__":
    main()

