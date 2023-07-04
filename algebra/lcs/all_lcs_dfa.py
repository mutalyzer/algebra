from collections import deque

from algebra.variants import Variant


class _Node:
    """Node class for (stretches of) matching symbols.
    FOR INTERNAL USE ONLY.
    """

    def __init__(self, row, col, length=0):
        self.row = row
        self.col = col
        self.length = length

        self.edges = []
        self.incoming = 0

    def __eq__(self, other):
        return (
            self.row == other.row
            and self.col == other.col
            and self.length == other.length
        )

    def __hash__(self):
        return hash((self.row, self.col, self.length))

    def __repr__(self):
        return f"{self.row, self.col, self.length}"


def to_dot(reference, source):
    """The LCS graph in Graphviz DOT format."""

    def traverse():
        # breadth-first traversal
        node_count = 0
        visited = {source: node_count}
        queue = deque([source])
        while queue:
            node = queue.popleft()

            if not node.edges:
                yield f'"s{visited[node]}"' + "[shape=doublecircle]"

            for succ, variant in node.edges:
                if succ not in visited:
                    node_count += 1
                    visited[succ] = node_count
                    queue.append(succ)
                yield (
                    f'"s{visited[node]}" -> "s{visited[succ]}"'
                    f' [label="{variant.to_hgvs(reference)}"];'
                )

    return (
        "digraph {\n"
        "    rankdir=LR\n"
        "    edge[fontname=monospace]\n"
        "    node[shape=circle]\n"
        "    si[shape=point]\n"
        "    si->" + f'"s{0}"' + "\n"
        "    " + "\n    ".join(traverse()) + "\n}"
    )


def lcs_graph_dfa(reference, observed, lcs_nodes):

    edges = []

    if not lcs_nodes or lcs_nodes == [[]]:
        source = _Node(0, 0, 0)
        if not reference and not observed:
            return source, edges
        sink = _Node(len(reference) + 1, len(observed) + 1, 0)
        variant = Variant(0, len(reference), observed)
        edges.append(variant)
        source.edges = [(sink, variant)]
        return source, edges

    if lcs_nodes[-1][-1].row + lcs_nodes[-1][-1].length == len(reference) + 1 and lcs_nodes[-1][-1].col + lcs_nodes[-1][-1].length == len(observed) + 1:
        sink = lcs_nodes[-1][-1]
        sink.length += 1
        for pred in lcs_nodes[-1][:-1]:
            if pred.row + pred.length - 1 < sink.row + sink.length - 1 and pred.col + pred.length - 1 < sink.col + sink.length - 1:
                variant = Variant(pred.row + pred.length - 1, sink.row + sink.length - 2, observed[pred.col + pred.length - 1:sink.col + sink.length - 2])
                pred.edges.append((sink, variant))
                edges.append(variant)
        sink.length -= 1
    else:
        sink = _Node(len(reference) + 1, len(observed) + 1, 1)
        for pred in lcs_nodes[-1]:
            if pred.row + pred.length - 1 < sink.row + sink.length - 1 and pred.col + pred.length - 1 < sink.col + sink.length - 1:
                variant = Variant(pred.row + pred.length - 1, sink.row + sink.length - 2, observed[pred.col + pred.length - 1:sink.col + sink.length - 2])
                pred.edges.append((sink, variant))
                edges.append(variant)

    lcs_pos = len(lcs_nodes) - 1

    while lcs_pos > 0:

        nodes = lcs_nodes[lcs_pos]

        idx_split = [False for e in lcs_nodes[lcs_pos - 1]]

        while nodes:
            node = nodes.pop(0)

            if not node.edges and node != sink:
                continue

            idx_pred = -1
            edge_added = False

            for i, pred in enumerate(lcs_nodes[lcs_pos - 1]):

                if pred.row + pred.length - 1 < node.row + node.length - 1 and pred.col + pred.length - 1 < node.col + node.length - 1:
                    variant = Variant(pred.row + pred.length - 1, node.row + node.length - 2, observed[pred.col + pred.length - 1:node.col + node.length - 2])
                    edges.append(variant)

                    if pred.incoming == lcs_pos:
                        idx_split[i] = False
                        upper_node = _Node(pred.row, pred.col, pred.length)
                        upper_node.edges = pred.edges + [(node, variant)]
                        pred.length += 1
                        lcs_nodes[lcs_pos - 1][i] = upper_node
                    else:
                        pred.edges.append((node, variant))

                    node.incoming = lcs_pos
                    edge_added = True
                    idx_pred = i

            if node.length > 1:
                idx_insert = idx_pred + 1
                node.length -= 1
                lcs_nodes[lcs_pos - 1].insert(idx_insert, node)
                idx_split.insert(idx_insert, edge_added)

        lcs_pos -= 1

    if lcs_nodes[0][0].row == 1 == lcs_nodes[0][0].col:
        source = lcs_nodes[0][0]
        source.row -= 1
        source.col -= 1
        successors = lcs_nodes[0][1:]
    else:
        source = _Node(0, 0, 1)
        successors = lcs_nodes[0]

    for succ in successors:
        if source.row + source.length - 1 < sink.row + sink.length - 1 and source.col + source.length - 1 < sink.col + sink.length - 1 and (succ == sink or succ.edges):
            variant = Variant(source.row + source.length - 1, succ.row + succ.length - 2, observed[source.col + source.length - 1:succ.col + succ.length - 2])
            source.edges.append((succ, variant))
            edges.append(variant)

    return source, edges


def traversal(root, atomics=False):
    """
    Traverse the LCS graph.
    """

    def traverse(node, path):
        if not node.edges:
            yield path

        for succ, variant in node.edges:
            if atomics:
                for atomic in variant.atomics():
                    yield from traverse(succ, path + atomic)
            else:
                yield from traverse(succ, path + [variant])

    yield from traverse(root, [])
