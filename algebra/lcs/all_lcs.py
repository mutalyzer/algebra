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
        self.pre_edges = []
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

    sink = _Node(len(reference) + 1, len(observed) + 1)
    source = _Node(0, 0)

    if reference == observed:
        source.edges = [(sink, [])]
        return source, []

    if not lcs_nodes or lcs_nodes == [[]]:
        variant = Variant(0, len(reference), observed)
        source.edges = [(sink, [variant])]
        return source, [variant]

    edges = []
    for node in lcs_nodes[-1]:
        offset = node.length
        variant = Variant(node.row + offset - 1, len(reference), observed[node.col + offset - 1:])
        if variant:
            node.edges = [(sink, [variant])]
            edges.append(variant)
        else:
            node.edges = [(sink, [])]

    for idx, nodes in enumerate(lcs_nodes[:0:-1]):
        lcs_pos = len(lcs_nodes) - idx - 1

        while nodes:
            node = nodes.pop(0)
            if not node.edges and not node.pre_edges:
                continue

            offset = node.length - 1
            for pred in nodes:
                if pred.length <= 1:
                    continue

                pred_offset = pred.length - 1
                if node.row + offset >= pred.row + pred_offset and node.col + offset >= pred.col + pred_offset:
                    if node.pre_edges:
                        split = _Node(node.row, node.col, node.length - 1)
                        split.edges = node.pre_edges + [(node, [])]
                        node.row += offset
                        node.col += offset
                        node.length = 1
                        node.pre_edges = []
                        offset = 0
                        lcs_nodes[lcs_pos - 1].append(split)

                    variant = Variant(pred.row + pred_offset - 1, node.row + offset - 1, observed[pred.col + pred_offset - 1:node.col + offset - 1])
                    pred.pre_edges.append((node, [variant]))
                    edges.append(variant)
                    node.incoming = lcs_pos

            for pred_idx, pred in enumerate(lcs_nodes[lcs_pos - 1]):
                pred_offset = pred.length
                if node.row + offset >= pred.row + pred_offset and node.col + offset >= pred.col + pred_offset:
                    if node.row + offset == pred.row + pred_offset and node.col + offset == pred.col + pred_offset:
                        continue

                    if pred.incoming == lcs_pos:
                        split = _Node(pred.row, pred.col, pred.length)
                        pred.row += pred_offset
                        pred.col += pred_offset
                        pred.length = 1
                        split.edges = [(pred, [])]
                        lcs_nodes[lcs_pos - 1][pred_idx] = split
                        pred = split
                    elif node.pre_edges:
                        split = _Node(node.row, node.col, node.length - 1)
                        split.edges = node.pre_edges + [(node, [])]
                        node.row += offset
                        node.col += offset
                        node.length = 1
                        node.pre_edges = []
                        offset = 0
                        lcs_nodes[lcs_pos - 1].append(split)

                    variant = Variant(pred.row + pred_offset - 1, node.row + offset - 1, observed[pred.col + pred_offset - 1:node.col + offset - 1])
                    pred.edges.append((node, [variant]))
                    edges.append(variant)
                    node.incoming = lcs_pos

            node.edges += node.pre_edges
            node.pre_edges = []

            if node.length > 1:
                node.length -= 1
                lcs_nodes[lcs_pos - 1].append(node)

    for node in lcs_nodes[0]:
        if node.edges:
            variant = Variant(0, node.row - 1, observed[:node.col - 1])
            if variant:
                source.edges.append((node, [variant]))
                edges.append(variant)
            else:
                source.edges.append((node, []))

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
        # depth-first traversal
        if not node.edges:
            yield path
            return

        for succ, variant in node.edges:
            if atomics and variant:
                for atomic in variant[0].atomics():
                    yield from traverse(succ, path + atomic)
            else:
                yield from traverse(succ, path + variant)

    yield from traverse(root, [])
