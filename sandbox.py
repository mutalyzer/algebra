import sys
from algebra.lcs import edit, lcs_graph, to_dot
from algebra.utils import random_sequence, random_variants
from algebra.variants import Parser, Variant, patch, to_hgvs


def reduce(reference, root):
    queue = [(root, None, None)]
    visited = set()
    dot = "digraph {\n"
    while queue:
        node, parent, in_variant = queue.pop(0)
        visited.add(node)
        print("pop", node)
        if not node.edges:
            dot += f'    "{parent.row}_{parent.col}" -> "{node.row}_{node.col}" [label="{to_hgvs(in_variant, reference, sequence_prefix=False)}"];\n'
            continue

        successors = {}
        for succ, out_variant in node.edges:
            if succ in visited:
                successors = {}
                continue

            length = len(out_variant[0]) if len(out_variant) else 0
            if length in successors:
                found = False
                for idx, value in enumerate(successors[length]):
                    if succ == value[0]:
                        # variant = Variant(min(value[1][0].start, out_variant[0].start), max(value[1][0].end, out_variant[0].end), value[1][0].sequence)
                        # successors[length][idx] = (succ, [variant])
                        print("repeat", variant)
                        found = True
                        break
                if not found:
                    successors[length].append((succ, out_variant))
            else:
                successors[length] = [(succ, out_variant)]

        for _, edges in sorted(successors.items(), reverse=True):
            for succ, variant in edges:
                queue.append((succ, node, variant))
                print("push", succ)

        if parent and successors:
            dot += f'    "{parent.row}_{parent.col}" -> "{node.row}_{node.col}" [label="{to_hgvs(in_variant, reference, sequence_prefix=False)}"];\n'

    return dot + "}"


def main():
    if len(sys.argv) < 3:
        reference = random_sequence(100, 100)
        variants = list(random_variants(reference, p=0.02, mu_deletion=1.5, mu_insertion=1.7))
    else:
        reference = sys.argv[1]
        variants = Parser(sys.argv[2]).hgvs()

    observed = patch(reference, variants)

    distance, lcs_nodes = edit(reference, observed)

    print(to_hgvs(variants, reference))
    print(distance)
    print(sum([len(level) for level in lcs_nodes]))

    root, _ = lcs_graph(reference, observed, lcs_nodes)

    print(to_dot(reference, root))
    print(reduce(reference, root))


if __name__ == "__main__":
    main()
