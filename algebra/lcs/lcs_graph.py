"""An efficient method of calculating all Longest Common Subsequence (LCS)
alignments. This method builds upon [1]_ by creating a compressed LCS
graph.

See Also
--------
algebra.lcs.edit_distance : Calculate only the simple edit distance.

References
----------
[1] S. Wu, U. Manber, G. Myers and W. Miller. "An O(NP) Sequence
Comparison Algorithm". In: Information Processing Letters, 35.6 (1990),
pp. 317-323.
"""


from collections import deque
from operator import attrgetter
from os.path import commonprefix
from ..variants import Variant, patch


class LCSgraph:
    """LCS graph class for storing all minimal alignments.

    Attributes
    ----------
    distance : int
        The simple edit distance.
    supremal : `Variant`
        The supremal variant representation.
    """

    class Node:
        """Node class for (stretches of) matching symbols."""
        def __init__(self, row, col, length):
            self.row = row
            self.col = col
            self.length = length

            self.edges = []

        def __eq__(self, other):
            return (self.row == other.row and self.col == other.col and
                    self.length == other.length)

        def __hash__(self):
            return hash((self.row, self.col, self.length))

        def __lt__(self, other):
            return (self.row, self.col) < (other.row, other.col)

        def __repr__(self):
            return f"{self.row, self.col}[{self.length}]"

    def __init__(self, reference, observed, shift=0, max_distance=None):
        self.distance, lcs_nodes = _lcs_nodes(reference, observed, shift, max_distance)
        self._source, self.supremal = _build_graph(reference, observed, lcs_nodes, shift)

    @classmethod
    def from_sequence(cls, reference, observed):
        """The supremal LCS graph for two sequences."""
        if reference == observed:
            return cls("", "")

        prefix_len, suffix_len = trim(reference, observed)
        return cls.from_variant(reference, [Variant(prefix_len, len(reference) - suffix_len, observed[prefix_len:len(observed) - suffix_len])])

    @classmethod
    def from_supremal(cls, reference, supremal):
        """The supremal LCS graph for a supremal variant."""
        return cls(reference[supremal.start:supremal.end], supremal.sequence, shift=supremal.start)

    @classmethod
    def from_variant(cls, reference, variants, offset=10):
        """Iteratively find the supremal LCS graph for an allele by
        repeatedly widening a range of influence.

        Parameters
        ----------
        reference : str
            The reference sequence.
        variants : list of `Variant`s
            The allele of interest.

        Other Parameters
        ----------------
        offset : int, optional
            The minimum offset around the variants.

        Returns
        -------
        graph : `LCSgraph`
            The LCS graph.
        """

        if not variants:
            return cls("", "")

        start = min(variants, key=attrgetter("start")).start
        end = max(variants, key=attrgetter("end")).end
        observed = patch(reference[start:end], [Variant(variant.start - start, variant.end - start, variant.sequence) for variant in variants])
        variant = Variant(start, end, observed)

        offset = max(offset, len(variant) // 2, 1)

        while True:
            start = max(0, variant.start - offset)
            end = min(len(reference), variant.end + offset)
            observed = reference[start:variant.start] + variant.sequence + reference[variant.end:end]

            graph = cls(reference[start:end], observed, shift=start)
            if graph.distance == 0:
                return graph

            if ((graph.supremal.start > start or graph.supremal.start == 0) and
                    (graph.supremal.end < end or graph.supremal.end == len(reference))):
                return graph

            offset *= 2

    def bfs_traversal(self, atomics=False):
        """Generate all (nodes and) edges in the LCS graph.

        Other Parameters
        ----------------
        atomics : bool, optional
            If set to `True` the variants are represented using separate
            deletions and insertions.

        Yields
        ------
        source node : `LCSgraph.Node`
            The source node.
        sink node : `LCSgraph.Node`
            The sink node.
        variant : list
            An edge, i.e., a sorted list of variants unique between the source
            and sink nodes.
        """

        visited = set()
        queue = deque([self._source])
        while queue:
            source = queue.popleft()
            if source in visited:
                continue

            for sink, variant in source.edges:
                if atomics:
                    for atomic in variant.atomics():
                        yield source, sink, atomic
                else:
                    yield source, sink, [variant]
                queue.append(sink)

            visited.add(source)

    def edges(self):
        """Set of all edges in the LCS graph."""
        return {edge[0] for *_, edge in self.bfs_traversal(atomics=False)}

    def nodes(self):
        """Generate all nodes in the LCS graph in topological order."""
        visited = {self._source}
        stack = [self._source]
        while stack:
            source = stack.pop()
            yield source

            for sink, _ in source.edges:
                visited.add(sink)
                stack.append(sink)

    def paths(self, atomics=False):
        """Traverse all paths (alignments) in the LCS graph.

        Parameters
        ----------
        atomics : bool, optional
            If set to `True` the variants are represented using separate
            deletions and insertions.

        Yields
        ------
        list
            A sorted list of variants representing a single LCS alignment.
        """

        def traversal(source, path):
            if not source.edges:
                yield path

            for child, variant in source.edges:
                if atomics:
                    for atomic in variant.atomics():
                        yield from traversal(child, path + atomic)
                else:
                    yield from traversal(child, path + [variant])

        return traversal(self._source, [])


def trim(lhs, rhs):
    """Find the lengths of the common prefix and common suffix between
    two sequences."""
    idx = len(commonprefix([lhs, rhs]))
    return idx, len(commonprefix([lhs[idx:][::-1], rhs[idx:][::-1]]))


def _lcs_nodes(reference, observed, shift=0, max_distance=None):
    """Calculate the simple edit distance between two strings (reference
    and observed) and construct a collection of LCS nodes.

    Other Parameters
    ----------------
    shift : int, optional
        Shift all variant positions with a given offset.
    max_distance : int, optional
        Stops the calculation if the simple edit distance exceeds this
        value.

    Raises
    ------
    ValueError
        If the calculation exceeds the optional maximum distance.

    Returns
    -------
    int
        The simple edit distance.
    list
        A collection of LCS nodes used to construct the LCS graph.

    See Also
    --------
    `_build_graph` : Constructs the LCS graph from LCS nodes.
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
                lcs_nodes[lcs_pos].append(LCSgraph.Node(match_row + shift, match_col + shift, row - match_row))
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
            lcs_nodes[lcs_pos].append(LCSgraph.Node(match_row + shift, match_col + shift, row - match_row))

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

        if max_distance and abs(delta) + 2 * (it - 1) > max_distance:
            raise ValueError("maximum distance exceeded")

    return abs(delta) + 2 * (it - 1), lcs_nodes[:max_lcs_pos + 1]


def _build_graph(reference, observed, lcs_nodes, shift=0):
    """Construct the compressed LCS graph from a collection of LCS nodes.
    `lcs_nodes` is destroyed during this process.

    Other Parameters
    ----------------
    shift : int, optional
        Shift all variant positions with a given offset.

    Returns
    -------
    `LCSgraph.Node`
        The source of the LCS graph.

    See Also
    --------
    `_lcs_nodes` : Calculates the LCS nodes.
    """

    if not lcs_nodes or lcs_nodes == [[]]:
        source = LCSgraph.Node(shift, shift, 0)
        if not reference and not observed:
            return source, Variant(0, 0, "")
        sink = LCSgraph.Node(len(reference), len(observed), 0)
        variant = Variant(shift, shift + len(reference), observed)
        source.edges = [(sink, variant)]
        return source, variant

    sink = lcs_nodes[-1][-1]
    if sink.row + sink.length == len(reference) + shift and sink.col + sink.length == len(observed) + shift:
        del lcs_nodes[-1][-1]
        sink.length += 1
    else:
        sink = LCSgraph.Node(len(reference) + shift, len(observed) + shift, 1)
    lcs_nodes.append([sink])

    max_sink = 0
    length = {}
    incoming = {}
    while len(lcs_nodes) > 1:
        while lcs_nodes[-1]:
            node = lcs_nodes[-1].pop(0)

            if node != sink and not node.edges:
                continue

            node_length = length.get(id(node), node.length)
            idx_parent = 0
            for idx, parent in enumerate(lcs_nodes[-2]):
                parent_length = length.get(id(parent), parent.length)

                if parent.row + parent_length < node.row + node_length and parent.col + parent_length < node.col + node_length:
                    variant = Variant(parent.row + parent_length, node.row + node_length - 1,
                                      observed[parent.col + parent_length - shift:node.col + node_length - 1 - shift])

                    if node == sink:
                        max_sink = max(max_sink, node.row + node_length - 1)

                    if incoming.get(id(parent), 0) == len(lcs_nodes):
                        split = LCSgraph.Node(parent.row, parent.col, parent.length)
                        split.edges = parent.edges + [(node, variant)]
                        lcs_nodes[-2][idx] = split
                        length[id(split)] = parent_length
                        parent.row += parent_length
                        parent.col += parent_length
                        parent.length -= parent_length
                    else:
                        parent.edges.append((node, variant))

                    incoming[id(node)] = len(lcs_nodes)
                    idx_parent = idx + 1

            if node_length > 1:
                lcs_nodes[-2].insert(idx_parent, node)
                length[id(node)] = node_length - 1

        del lcs_nodes[-1]

    source = lcs_nodes[0][0]
    if lcs_nodes[0][0].row == lcs_nodes[0][0].col == shift:
        del lcs_nodes[0][0]
    else:
        source = LCSgraph.Node(shift, shift, 0)

    for node in lcs_nodes[0]:
        if node != sink and not node.edges:
            continue

        node_length = length.get(id(node), node.length)
        if source.row < node.row + node_length and source.col < node.col + node_length:
            variant = Variant(source.row, node.row + node_length - 1,
                              observed[source.col - shift:node.col + node_length - 1 - shift])

            if node == sink:
                max_sink = max(max_sink, node.row + node_length - 1)

            source.edges.append((node, variant))

    source_offset = min((edge.start for _, edge in source.edges), default=shift) - shift
    source.row += source_offset
    source.col += source_offset
    source.length -= source_offset
    sink.length -= sink.row + sink.length - max_sink

    return source, Variant(source.row, sink.row + sink.length, observed[source.col - shift:sink.col + sink.length - shift])
