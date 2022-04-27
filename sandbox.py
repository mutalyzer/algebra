from itertools import product
import sys
from algebra.lcs.wupp import edit, lcs_graph
from algebra.utils import random_sequence
from algebra.variants.variant import Variant


def collapse(variants):
    collapsed = []
    deletion = Variant(0, 0)
    insertion = None
    for variant in sorted(variants, key=lambda v: (v.start, v.end, len(v.sequence))):
        if variant.start < variant.end and len(variant.sequence) == 0:
            if variant.start <= deletion.end:
                deletion.end = variant.end
            else:
                if deletion:
                    collapsed.append(deletion)
                deletion = variant
        else:
            if insertion is not None and variant.start >= insertion.start and variant.start < insertion.end:
                insertion.sequence = "".join(set(insertion.sequence) | set(variant.sequence))
                insertion.end = max(insertion.end, variant.end)
            else:
                if insertion:
                    collapsed.append(insertion)
                insertion = Variant(variant.start, variant.end, "".join(set(variant.sequence)))
    if deletion:
        collapsed.append(deletion)
    if insertion:
        collapsed.append(insertion)

    return collapsed


def ops_set(edges):
    def explode(variant):
        for pos in range(variant.start, variant.end):
            yield Variant(pos, pos + 1)
        for pos in range(variant.start, variant.end + 1):
            for symbol in variant.sequence:
                yield Variant(pos, pos, symbol)

    ops = set()
    for edge in edges:
        ops.update(explode(edge))

    return ops


def check(reference, observed):
    _, lcs_nodes = edit(reference, observed)
    _, edges = lcs_graph(reference, observed, lcs_nodes)

    print("EDGES:", [variant.to_hgvs(reference) for variant in sorted(edges, key=lambda v: (v.start, v.end, len(v.sequence)))])

    collapsed_edges = collapse(edges)
    print("COLL :", [variant.to_hgvs(reference) for variant in sorted(collapsed_edges, key=lambda v: (v.start, v.end, len(v.sequence)))])

    ops = ops_set(edges)
    collapsed_ops = ops_set(collapsed_edges)

    print("OPS  :", [variant.to_hgvs(reference) for variant in sorted(ops, key=lambda v: (v.start, v.end, len(v.sequence)))])

    print(ops - collapsed_ops)
    print(collapsed_ops - ops)

    assert ops == collapsed_ops


def main():
    if len(sys.argv) < 3:
        reference = random_sequence(15)
        observed = random_sequence(15)
    else:
        reference = sys.argv[1]
        observed = sys.argv[2]

    print(f'"{reference}" "{observed}"')
    check(reference, observed)


if __name__ == "__main__":
    main()
