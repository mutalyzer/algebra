"""An efficient method of calculating all Longest Common Subsequence (LCS)
alignments. This method builds upon [1]_ by creating a compressed LCS
graph in two phases.

Phase 1 (`edit`) calculates the simple edit distance and a collection of
LCS nodes.
Phase 2 (`lcs_graph`) creates a directed acyclic (multi) graph of the LCS
nodes. Every path in this graph is a unique and distinct embedding of an
LCS and consequently a unique variant representation.

See Also
--------
algebra.lcs.distance_only : Calculates only the simple edit distance.

References
----------
[1] S. Wu, U. Manber, G. Myers and W. Miller. "An O(NP) Sequence
Comparison Algorithm". In: Information Processing Letters, 35.6 (1990),
pp. 317-323.
"""


from ..variants import Variant


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
        return (self.row == other.row and self.col == other.col and
                self.length == other.length)

    def __hash__(self):
        return hash((self.row, self.col, self.length))

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
                lcs_pos = ((row + col) - (abs(delta) + 2 * it - abs((len(reference) - row) - (len(observed) - col)))) // 2 - 1
                max_lcs_pos = max(lcs_pos, max_lcs_pos)
                lcs_nodes[lcs_pos].append(_Node(match_row + 1, match_col + 1, row - match_row))
                matching = False
            row += 1
            col += 1

        steps = end + 1
        if not matching:
            match_row = row
            match_col = col
        while row < len(reference) and col < len(observed) and reference[row] == observed[col]:
            matching = True
            row += 1
            col += 1
            steps += 1
        if matching:
            lcs_pos = ((row + col) - (abs(delta) + 2 * it - abs((len(reference) - row) - (len(observed) - col)))) // 2 - 1
            max_lcs_pos = max(lcs_pos, max_lcs_pos)
            lcs_nodes[lcs_pos].append(_Node(match_row + 1, match_col + 1, row - match_row))

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

        for idx in range(upper + it, delta, - 1):
            diagonals[offset + idx] = expand(idx)

        diagonals[offset + delta] = expand(delta)
        it += 1

    return abs(delta) + 2 * (it - 1), lcs_nodes[:max_lcs_pos + 1]


def lcs_graph(reference, observed, lcs_nodes):
    """Construct the compressed LCS graph. `lcs_nodes` is destroyed
    during this process.

    Returns
    -------
    `_Node` (opaque data type)
        The root of the LCS graph.
    list
        A list of edges (variants) in the LCS graph.

    See Also
    --------
    `edit` : Calculates the LCS nodes.
    `traversal` : Traverses the LS graph.
    `algebra.utils.to_dot` : Graphviz DOT format.
    """


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

        idx_split = [False for _ in lcs_nodes[lcs_pos - 1]]

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
    """Traverse the LCS graph.

    Parameters
    ----------
    atomics : bool, optional
        If set to `True` the variants are represented using separate
        deletions and insertions.

    Yields
    ------
    list
        A sorted list of variants representing an LCS alignment.

    See Also
    --------
    `lcs_graph` : Constructs the LCS graph.
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
