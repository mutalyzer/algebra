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
        self.incoming = -2

    def __eq__(self, other):
        return (
            self.row == other.row
            and self.col == other.col
            and self.length == other.length
        )

    def __hash__(self):
        return id(self)

    def __repr__(self):
        return f"{self.row, self.col, self.length}"


def edit(reference, observed):
    """Calculate the simple edit distance between two strings and
    construct a collection of LCS nodes.

    Returns
    -------
    int
        The simple edit distance.
    list
        A collection of LCS nodes in order to construct the LCS graph.

    See Also
    --------
    `lcs_graph` : Constructs the LCS graph from LCS nodes.
    """

    def expand(idx):
        nonlocal max_lcs_pos
        start = diagonals[offset + idx]
        if idx > 0:
            row = start
            col = row + idx
            end = max(diagonals[offset + idx - 1] - 1, diagonals[offset + idx + 1])
        elif idx < 0:
            col = start
            row = col - idx
            end = max(diagonals[offset + idx - 1], diagonals[offset + idx + 1] - 1)
        else:
            row = start
            col = row + idx
            end = max(diagonals[offset + idx - 1], diagonals[offset + idx + 1])

        matching = False
        match_row = 0
        match_col = 0
        for _ in range(start, end):
            if reference[row] == observed[col]:
                if not matching:
                    match_row = row
                    match_col = col
                matching = True
            elif matching:
                lcs_pos = (
                    (row + col)
                    - (
                        abs(delta)
                        + 2 * it
                        - abs((len(reference) - row) - (len(observed) - col))
                    )
                ) // 2 - 1
                max_lcs_pos = max(lcs_pos, max_lcs_pos)
                lcs_nodes[lcs_pos].append(
                    _Node(match_row + 1, match_col + 1, row - match_row)
                )
                matching = False
            row += 1
            col += 1

        steps = end + 1
        if not matching:
            match_row = row
            match_col = col
        while (
            row < len(reference)
            and col < len(observed)
            and reference[row] == observed[col]
        ):
            matching = True
            row += 1
            col += 1
            steps += 1
        if matching:
            lcs_pos = (
                (row + col)
                - (
                    abs(delta)
                    + 2 * it
                    - abs((len(reference) - row) - (len(observed) - col))
                )
            ) // 2 - 1
            max_lcs_pos = max(lcs_pos, max_lcs_pos)
            lcs_nodes[lcs_pos].append(
                _Node(match_row + 1, match_col + 1, row - match_row)
            )

        return steps

    lcs_nodes = [[] for _ in range(min(len(reference), len(observed)))]
    max_lcs_pos = 0

    delta = len(observed) - len(reference)
    offset = len(reference) + 1
    diagonals = [0] * (len(reference) + len(observed) + 3)
    it = 0

    if delta >= 0:
        lower = 0
        upper = delta
    else:
        lower = delta
        upper = 0

    while diagonals[offset + delta] <= max(len(reference), len(observed)) - abs(delta):
        for idx in range(lower - it, delta):
            diagonals[offset + idx] = expand(idx)

        for idx in range(upper + it, delta, -1):
            diagonals[offset + idx] = expand(idx)

        diagonals[offset + delta] = expand(delta)
        it += 1

    return abs(delta) + 2 * (it - 1), lcs_nodes[: max_lcs_pos + 1]


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
                yield (f'"s{visited[node]}" -> "s{visited[succ]}"'
                       f' [label="{variant.to_hgvs(reference)}"];')

    return (
        "digraph {\n"
        "    rankdir=LR\n"
        "    edge[fontname=monospace]\n"
        "    node[shape=circle]\n"
        "    si[shape=point]\n"
        "    si->" + f'"s{0}"' + "\n"
        "    " + "\n    ".join(traverse()) + "\n}"
    )


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
    start = pred.row + pred.length - 1
    end = succ.row + succ.length - 2
    seq_start = pred.col + pred.length - 1
    seq_end = succ.col + succ.length - 2
    return Variant(start, end, observed[seq_start:seq_end])


def variant_possible(pred, succ):
    return (
        pred.row + pred.length - 1 < succ.row + succ.length - 1
        and pred.col + pred.length - 1 < succ.col + succ.length - 1
    )


def should_merge_sink(lcs_nodes, sink):
    last_node = lcs_nodes[-1][-1]
    return (
        last_node.row + last_node.length == sink.row
        and lcs_nodes[-1][-1].col + last_node.length == sink.col
    )


def get_node_with_length(node):
    return node.row + node.length - 1, node.col + node.length - 1


def lcs_graph_mdfa(reference, observed, lcs_nodes):
    sink = _Node(len(reference) + 1, len(observed) + 1, 1)
    source = _Node(0, 0, 1)
    print_lcs_nodes(lcs_nodes)

    edges = []

    if not reference and not observed:
        return None, []

    if not lcs_nodes or lcs_nodes == [[]]:
        variant = Variant(0, len(reference), observed)
        edges.append(variant)
        source.edges = [(sink, variant)]
        return source, edges

    predecessors = lcs_nodes[-1]

    if should_merge_sink(lcs_nodes, sink):
        print("should merge sink")
        lcs_nodes[-1][-1].length += 1
        sink = lcs_nodes[-1][-1]
        predecessors = lcs_nodes[-1][:-1]

    successors = [sink]
    lcs_pos = len(lcs_nodes) - 1


    while lcs_pos >= -1:
        print(f"\nfor level {lcs_pos} to {lcs_pos + 1}")
        print("----------------")
        print([(n, n.incoming) for n in predecessors])
        print([(n, n.incoming) for n in successors])
        print("----------------")

        while successors:
            node = successors.pop(0)
            print(f"\n- working with node ({node}, {node.incoming}) as {get_node_with_length(node)}")

            if not node.edges and node != sink:
                print(f"- skipped {node}")
                continue

            idx_pred = -1
            for i, pred in enumerate(predecessors):
                print(f" - predecessor node ({pred}, {pred.incoming}) as {get_node_with_length(pred)}")
                if variant_possible(pred, node):
                    variant = get_variant(pred, node, observed)

                    edges.append(variant)

                    if pred.incoming == lcs_pos:
                        print(f"\n  - inversion")
                        upper_node = _Node(pred.row, pred.col, pred.length)
                        upper_node.edges = pred.edges + [(node, variant)]
                        predecessors[i] = upper_node
                        predecessor = upper_node
                        print(f"   - split predecessor")
                        print(f"   - upper node: ({upper_node}, {upper_node.incoming}); edges: {upper_node.edges}\n")
                    else:
                        predecessor = pred
                        predecessor.edges.append((node, variant))
                    node.incoming = lcs_pos
                    idx_pred = i

                    print(f"  - added edge: {variant.to_hgvs(reference)}")
                    print(f"  - predecessor edges: {[((e[0], e[0].incoming), e[1].to_hgvs(reference)) for e in predecessor.edges]}")
                    print(f"  - set incoming for working node to: {node.incoming}")

            print(" - check if current node should be added to the predecessors")
            if node.length > 1:
                node.length -= 1
                predecessors.insert(idx_pred+1, node)
                print(f"  - added ({node}, {node.incoming}) at index {idx_pred + 1} in the predecessors:\n    {[(n, n.incoming) for n in predecessors]}")

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

    return source, edges


def traversal(root):
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
            if not variant:
                raise ValueError
            yield from traverse(succ, path + [variant])

    yield from traverse(root, [])
