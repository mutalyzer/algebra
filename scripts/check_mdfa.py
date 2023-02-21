from algebra.lcs import edit, lcs_graph, traversal
from algebra.lcs.all_lcs_mdfa import edit as edit_mdfa
from algebra.lcs.all_lcs_mdfa import lcs_graph_mdfa
from algebra.utils import random_sequence, random_variants
from algebra.variants import patch
from algebra.variants.variant import to_hgvs


def traversal_no_variant(root):
    """
    Traverse the LCS graph.
    """
    def traverse(node, path):
        if not node.edges:
            yield path

        variants = [str(e[1]) for e in node.edges]
        if len(variants) != len(set(variants)):
            raise ValueError
        for succ, variant in node.edges:
            empty_variants = [v for v in variant if len(v) == 0]
            if empty_variants != [] or variant == []:
                raise ValueError
            yield from traverse(succ, path + variant)

    yield from traverse(root, [])


def get_minimal_mdfa(reference, observed):
    _, lcs_nodes_mdfa = edit_mdfa(reference, observed)
    source_mdfa, sink_mdfa = lcs_graph_mdfa(reference, observed, lcs_nodes_mdfa)
    return traversal_no_variant(source_mdfa)


def get_minimal(reference, observed):
    _, lcs_nodes = edit(reference, observed)
    source, sink = lcs_graph(reference, observed, lcs_nodes)

    return traversal(source)


def compare_minimal(reference, observed):
    minimal_mdfa = sorted(
        [to_hgvs(vars, reference) for vars in get_minimal_mdfa(reference, observed)],
        key=str,
    )
    minimal = sorted(
        [to_hgvs(vars, reference) for vars in get_minimal(reference, observed)], key=str
    )

    print(minimal_mdfa)
    print(minimal)

    return minimal_mdfa == minimal


def main():

    equal = True
    count = 0
    while equal:
        print("\n\n=====================")
        reference = random_sequence(30, 1)
        observed = random_sequence(30, 1)
        # observed = patch(reference, list(random_variants(reference)))
        print(reference, observed)
        print(observed)
        print("===========")
        equal = compare_minimal(reference, observed)
        print("===========")
        print(equal)
        print("=====================\n\n")
        count += 1
        print("count:", count)

        if equal:
            # break
            continue
        else:
            print("error")
            break


if __name__ == "__main__":
    main()
