from algebra.lcs import edit, lcs_graph, traversal
from algebra.lcs.all_lcs_mdfa import edit as edit_mdfa
from algebra.lcs.all_lcs_mdfa import lcs_graph_mdfa
from algebra.utils import random_sequence, random_variants
from algebra.variants import patch
from algebra.variants.variant import to_hgvs


def get_minimal_mdfa(reference, observed):
    _, lcs_nodes_mdfa = edit_mdfa(reference, observed)
    source_mdfa, sink_mdfa = lcs_graph_mdfa(reference, observed, lcs_nodes_mdfa)

    return traversal(source_mdfa)


def get_minimal(reference, observed):
    _, lcs_nodes = edit(reference, observed)
    source, sink = lcs_graph(reference, observed, lcs_nodes)

    return traversal(source)


def extract(reference, observed):
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
    while equal:
        print("\n\n=====================")
        reference = random_sequence(7, 4)
        observed = patch(reference, list(random_variants(reference)))
        equal = extract(reference, observed)
        print(reference, observed)
        print(observed)
        print(equal)
        print("=====================")
        if equal:
            continue
        else:
            print("error")
            break


if __name__ == "__main__":
    main()
