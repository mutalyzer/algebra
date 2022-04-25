from ..variants.variant import Variant, to_hgvs


class Node:
    def __init__(self, row, col, length=0):
        self.row = row
        self.col = col
        self.length = length

        self.edges = []
        self.pre_edges = []
        self.incoming = 0

    def __repr__(self):
        return f"{self.row, self.col, self.length}"


def edit(reference, observed):
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

        active = False
        match_row = 0
        match_col = 0
        for _ in range(start, end):
            if reference[row] == observed[col]:
                if not active:
                    match_row = row
                    match_col = col
                active = True
            elif active:
                lcs_pos = ((row + col) - (abs(delta) + 2 * it - abs((len(reference) - row) - (len(observed) - col)))) // 2 - 1
                max_lcs_pos = max(lcs_pos, max_lcs_pos)
                lcs_nodes[lcs_pos].append(Node(match_row + 1, match_col + 1, row - match_row))
                active = False
            row += 1
            col += 1

        steps = end + 1
        if not active:
            match_row = row
            match_col = col
        while row < len(reference) and col < len(observed) and reference[row] == observed[col]:
            active = True
            row += 1
            col += 1
            steps += 1
        if active:
            lcs_pos = ((row + col) - (abs(delta) + 2 * it - abs((len(reference) - row) - (len(observed) - col)))) // 2 - 1
            max_lcs_pos = max(lcs_pos, max_lcs_pos)
            lcs_nodes[lcs_pos].append(Node(match_row + 1, match_col + 1, row - match_row))

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
    sink = Node(len(reference) + 1, len(observed) + 1)

    source = Node(0, 0)
    if lcs_nodes == [[]]:
        source.edges = [(sink, [Variant(0, len(reference), observed)])]
        return source

    for node in lcs_nodes[-1]:
        offset = node.length
        variant = Variant(node.row + offset - 1, len(reference), observed[node.col + offset - 1:])
        node.edges = [(sink, [variant] if variant else [])]

    for idx, nodes in enumerate(lcs_nodes[:0:-1]):
        lcs_pos = len(lcs_nodes) - idx - 1

        while nodes:
            node = nodes.pop(0)
            if not node.edges and not node.pre_edges:
                continue

            offset = node.length - 1
            for pred_idx, pred in enumerate(lcs_nodes[lcs_pos - 1]):
                pred_offset = pred.length
                if node.row + offset >= pred.row + pred_offset and node.col + offset >= pred.col + pred_offset:
                    if pred.incoming == lcs_pos:
                        split = Node(pred.row, pred.col, pred.length)
                        pred.row += pred_offset
                        pred.col += pred_offset
                        pred.length = 1
                        split.edges = [(pred, [])]
                        lcs_nodes[lcs_pos - 1][pred_idx] = split
                        pred = split
                    elif node.pre_edges:
                        split = Node(node.row, node.col, node.length - 1)
                        split.edges = node.pre_edges + [(node, [])]
                        node.row += offset
                        node.col += offset
                        node.length = 1
                        node.pre_edges = []
                        offset = 0
                        lcs_nodes[lcs_pos - 1].append(split)

                    variant = Variant(pred.row + pred_offset - 1, node.row + offset - 1, observed[pred.col + pred_offset - 1:node.col + offset - 1])
                    pred.edges.append((node, [variant]))
                    node.incoming = lcs_pos

            for pred in nodes:
                if pred.length <= 1:
                    continue

                pred_offset = pred.length - 1
                if node.row + offset >= pred.row + pred_offset and node.col + offset >= pred.col + pred_offset:
                    if node.pre_edges:
                        split = Node(node.row, node.col, node.length - 1)
                        split.edges = node.pre_edges + [(node, [])]
                        node.row += offset
                        node.col += offset
                        node.length = 1
                        node.pre_edges = []
                        offset = 0
                        lcs_nodes[lcs_pos - 1].append(split)

                    variant = Variant(pred.row + pred_offset - 1, node.row + offset - 1, observed[pred.col + pred_offset - 1:node.col + offset - 1])
                    pred.pre_edges.append((node, [variant]))
                    node.incoming = lcs_pos

            node.edges += node.pre_edges
            node.pre_edges = []

            if node.length > 1:
                node.length -= 1
                lcs_nodes[lcs_pos - 1].append(node)

    for node in lcs_nodes[0]:
        if node.edges:
            variant = Variant(0, node.row - 1, observed[:node.col - 1])
            source.edges.append((node, [variant] if variant else []))

    return source


def traversal(root, atomics=False):
    def traverse(node, path):
        if not node.edges:
            yield path
            return

        for succ, variant in node.edges:
            if atomics and len(variant) > 0:
                for atomic in variant[0].atomics():
                    yield from traverse(succ, path + atomic)
            else:
                yield from traverse(succ, path + variant)

    yield from traverse(root, [])


def to_dot(reference, root):
    def nodes_and_edges():
        queue = [root]
        visited = [root]
        while queue:
            node = queue.pop(0)
            yield f'"{node.row}_{node.col}" [label="{node.row, node.col}"];'
            for succ, variant in node.edges:
                yield f'"{node.row}_{node.col}" -> "{succ.row}_{succ.col}" [label="{to_hgvs(variant, reference, sequence_prefix=False)}"];'
                if succ not in visited:
                    visited.append(succ)
                    queue.append(succ)

    return "digraph {\n    " + "\n    ".join(nodes_and_edges()) + "\n}"
