from algebra.lcs.efficient import edit as edit_old, build, traversal as traversal_old
from algebra.lcs.wupp import edit as edit_new, lcs_graph, traversal as traversal_new
from algebra.variants.variant import Variant, to_hgvs
import sys


def explore(reference, observed, root):
    def explode(variant):
        for pos in range(variant.start, variant.end):
            yield Variant(pos, pos + 1)
        for pos in range(variant.start, variant.end + 1):
            for symbol in variant.sequence:
                yield Variant(pos, pos, symbol)

    ops = set()
    queue = [root]
    visited = [root]
    while queue:
        node = queue.pop(0)
        max_succ = (0, 0)
        for succ, variant in node.edges:
            if succ not in visited:
                queue.append(succ)
                visited.append(succ)
            max_succ = (max(max_succ[0], succ.row), max(max_succ[1], succ.col))
        if max_succ != (0, 0):
            variant = Variant(node.row, max_succ[0] - 1, observed[node.col:max_succ[1] - 1])
            if variant:
                print(variant.to_hgvs(reference), {var.to_hgvs() for var in explode(variant)})
                ops |= {var for var in explode(variant)}
    return ops


def main():
    if len(sys.argv) < 3:
        reference = "CATATATCG"
        observed = "CTTATAGCAT"
    else:
        reference = sys.argv[1]
        observed = sys.argv[2]

    distance_old, lcs_nodes_old = edit_old(reference, observed)
    ops_old, graph_old = build(lcs_nodes_old, reference, observed)

    hgvs = set()
    for op in ops_old:
        if op[1] == "ins":
            variant = Variant(op[0], op[0], op[2])
        else:
            variant = Variant(op[0] - 1, op[0])
        hgvs.add(variant.to_hgvs())

    print(sorted(hgvs))
    print(len(hgvs))
    #for path in traversal_old(reference, observed, graph_old, atomics=True):
    #    print(to_hgvs(path, reference, sequence_prefix=False))

    distance_new, lcs_nodes_new = edit_new(reference, observed)
    graph_new = lcs_graph(reference, observed, lcs_nodes_new)

    ops_new = explore(reference, observed, graph_new)
    print({op.to_hgvs() for op in ops_new})
    print(len(ops_new))


if __name__ == "__main__":
    main()
