from argparse import ArgumentParser

from algebra.lcs import edit
from algebra.lcs.all_lcs import _Node
from algebra.utils import to_dot
from algebra.variants import Variant


def print_lcs_nodes(lcs_nodes):
    print("\n=========")
    print("lcs nodes")
    print("=========")
    for idx, level in enumerate(lcs_nodes):
        print(f"level: {idx}")
        print("--------")
        for node in level:
            print(node, node.edges)
        print("--------")
    print("=========\n")


def get_variant(pred, succ, observed):
    return Variant(pred.row, succ.row - 1, observed[pred.col : succ.col - 1])


def get_variant_offset(pred, succ, observed):
    start = pred.row + pred.length - 1
    end = succ.row + succ.length - 2
    seq_start = pred.col + pred.length - 1
    seq_end = succ.col + succ.length - 2
    return Variant(start, end, observed[seq_start:seq_end])


def variant_possible_offset(pred, succ):
    return (
        pred.row + pred.length - 1 < succ.row + succ.length - 1
        and pred.col + pred.length - 1 < succ.col + succ.length - 1
    )


def get_edge(pred, succ, observed):
    return succ, get_variant(pred, succ, observed)


def variant_possible(pred, succ):
    return pred.row < succ.row and pred.col < succ.col


def print_edge(pred, succ, reference, observed):
    if variant_possible(pred, succ):
        edge = get_edge(pred, succ, observed)[1]
        print(f"{pred} -> {succ}: {edge} -> {edge.to_hgvs(reference)}")
    else:
        print(f"{pred} -> {succ}: no edge possible")


def should_merge_sink(lcs_nodes, sink):
    return (
        lcs_nodes[-1][-1].row + 1 == sink.row and lcs_nodes[-1][-1].col + 1 == sink.col
    )


def get_node_with_length(node):
    return node.row + node.length - 1, node.col + node.length - 1


def split_node(node):
    lower_node = _Node(node.row + node.length, node.col + node.length)
    lower_node.edges = node.edges
    lower_node.pre_edges = node.pre_edges
    node.edges = [(lower_node, [])]
    node.pre_edges = []
    for pre_edge in lower_node.pre_edges:
        node_to_update = pre_edge[0]
        for i, e in enumerate(node_to_update.edges):
            if e[0] == node:
                node_to_update.edges[i] = (lower_node, e[1])


def lcs_graph_mdfa(reference, observed, lcs_nodes):
    sink = _Node(len(reference) + 1, len(observed) + 1, 1)
    source = _Node(0, 0, 1)

    print_lcs_nodes(lcs_nodes)

    print_edge(_Node(3, 3), _Node(3, 6), reference, observed)
    print_edge(_Node(3, 3), _Node(4, 7), reference, observed)

    predecessors = lcs_nodes[-1]

    if should_merge_sink(lcs_nodes, sink):
        lcs_nodes[-1][-1].length += 1
        sink = lcs_nodes[-1][-1]
        predecessors = lcs_nodes[-1][:-1]

    successors = [sink]

    lcs_pos = len(lcs_nodes) - 1

    print("predecessors")
    print(predecessors)
    print("successors")
    print(successors)

    while predecessors:
        print(f"\nfor level {lcs_pos} to {lcs_pos + 1}")
        print("----------------")
        print(predecessors)
        print(successors)
        print("----------------")

        while successors:
            node = successors.pop(0)
            print(f"\n- working with node {node} as {get_node_with_length(node)}")
            if not node.edges and node != sink:
                continue

            for pred in predecessors:
                print(f" - predecessor node {pred} as {get_node_with_length(pred)}")
                if variant_possible_offset(pred, node):
                    print(f"  - node edges: {node.edges}")
                    print(f"  - node pre edges: {node.pre_edges}")
                    print(f"  - pred node edges: {pred.edges}")
                    print(f"  - pred node pre edges: {pred.pre_edges}")
                    variant = get_variant_offset(pred, node, observed)
                    for pre_edge in pred.pre_edges:
                        if pre_edge[1][0].end == variant.end:
                            print(f"\n!!!possible inversion!!!\nshould split pred {pred}")
                            split_node(pred)
                    pred.edges.append((node, [variant]))
                    node.pre_edges.append((pred, [variant]))
                    print(f"  - added {variant}, {variant.to_hgvs(reference)}")

            print(" - check if current node should be added to the predecessors")
            if node.length > 1:
                node.length -= 1
                predecessors = sorted(predecessors + [node], key=lambda n: (n.row, n.col))
                print(f"  - added with new predecessors {predecessors}")


        lcs_pos -= 1
        successors = predecessors

        if lcs_pos > -1:
            predecessors = lcs_nodes[lcs_pos]
        elif lcs_pos == -1:
            if successors[0].row == 1 and successors[0].col == 1:
                source = successors[0]
                source.row = 0
                source.col = 0
                successors.pop(0)
            predecessors = [source]
        else:
            break

    return source, None


def extract(reference, observed):
    _, lcs_nodes = edit(reference, observed)
    root, _ = lcs_graph_mdfa(reference, observed, lcs_nodes)

    return to_dot(reference, root)


def main():
    # print(extract("TTT", "TTTTAT"))
    parser = ArgumentParser(description="Extract variants from line-based input")
    parser.add_argument("ref", type=str, help="reference sequence")
    parser.add_argument("obs", type=str, help="observed sequence")

    args = parser.parse_args()

    print(extract(args.ref, args.obs))


if __name__ == "__main__":
    main()