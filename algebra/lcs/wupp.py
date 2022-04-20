from ..variants.variant import Variant, to_hgvs


class Node:
    def __init__(self, row, col, len=0):
        self.row = row
        self.col = col
        self.len = len

        self.edges = []

        # TODO: init?!
        self.incoming = 0
        self.outgoing = 0

        self.internal = []

    def __repr__(self):
        # return f"{self.row, self.col, self.len, self.incoming, [(a, to_hgvs(b)) for a,b in self.edges]}"
        return f"{self.row, self.col, self.len, self.incoming, len(self.edges)}"


def edit(reference, observed):
    def lcs_idx(row, col):
        return ((row + col) - (abs(delta) + 2 * it - abs((len(reference) - row) - (len(observed) - col)))) // 2 - 1

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

        # print(f"    start at {row, col}")
        active = False
        match_row = 0
        match_col = 0
        for _ in range(start, end):
            # print(f"        implicit {row, col}")
            matrix[row + 1][col + 1] = abs(delta) + 2 * it
            if reference[row] == observed[col]:
                # print(f"{row + 1, col + 1} implicit match: {reference[row]}")
                if not active:
                    match_row = row
                    match_col = col
                active = True
            elif active:
                print(f"{match_row + 1, match_col + 1}, {row + 1, col + 1} match")
                lcs_pos = lcs_idx(row, col)
                max_lcs_pos = max(lcs_pos, max_lcs_pos)
                print('lcs pos', lcs_pos)
                lcs_nodes[lcs_pos].append(Node(match_row + 1, match_col + 1, row - match_row))
                active = False
            row += 1
            col += 1

        # print(f"    continue with {row, col}")
        steps = end + 1
        if not active:
            match_row = row
            match_col = col
        while row < len(reference) and col < len(observed) and reference[row] == observed[col]:
            active = True
            # print(f"        match {row, col}")
            matrix[row + 1][col + 1] = abs(delta) + 2 * it
            # print(f"{row + 1, col + 1} match: {reference[row]}")
            row += 1
            col += 1
            steps += 1
        if active:
            print(f"{match_row + 1, match_col + 1}, {row + 1, col + 1} match")
            lcs_pos = lcs_idx(row, col)
            max_lcs_pos = max(lcs_pos, max_lcs_pos)
            print('lcs pos', lcs_pos)
            lcs_nodes[lcs_pos].append(Node(match_row + 1, match_col + 1, row - match_row))

        # print(f"    last checked {row, col}")
        matrix[row + 1][col + 1] = abs(delta) + 2 * it + 2
        return steps

    matrix = [[None for _ in range(len(observed) + 2)] for _ in range(len(reference) + 2)]
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
            # print(f"extend diagonal {idx} (lower)")
            diagonals[offset + idx] = expand(idx)
        for idx in range(upper + it, delta, - 1):
            # print(f"extend diagonal {idx} (upper)")
            diagonals[offset + idx] = expand(idx)
        # print(f"extend diagonal {delta} (delta)")
        diagonals[offset + delta] = expand(delta)

        # for row in matrix:
        #     print(row)
        # print(diagonals)

        it += 1

    return abs(delta) + 2 * (it - 1), [r[:-1] for r in matrix[:-1]], lcs_nodes[:max_lcs_pos + 1]


def lcs_graph(reference, observed, lcs_nodes):
    sink = Node(len(reference) + 1, len(observed) + 1)

    for level in lcs_nodes:
        print('level', level)

    source = Node(0, 0)
    if lcs_nodes == [[]]:
        source.edges = [(sink, [Variant(0, len(reference), observed)])]
        return source

    for node in lcs_nodes[-1]:
        offset = node.len - 1
        variant = Variant(node.row + offset, len(reference), observed[node.col + offset:])
        node.edges = [(sink, [variant] if variant else [])]

    for idx, nodes in enumerate(lcs_nodes[:0:-1]):
        lcs_pos = len(lcs_nodes) - idx - 1

        while len(nodes) > 0:
            node = nodes.pop(0)

            if not node.edges and not node.internal:
                continue

            offset = node.len - 1

            for pred in nodes:
                if pred.len <= 1:
                    continue

                pred_offset = pred.len - 2

                if node.row + offset > pred.row + pred_offset and node.col + offset > pred.col + pred_offset:
                    print(f"S{node.row + offset, node.col + offset} vs {pred.row + pred_offset, pred.col + pred_offset}")

                    variant = Variant(pred.row + pred_offset, node.row + offset - 1, observed[pred.col + pred_offset:node.col + offset - 1])
                    print(variant.to_hgvs(reference))

                    if pred.incoming == lcs_pos:
                        print(f"Split incoming {pred}")

                        split = Node(pred.row, pred.col, pred.len - 1)

                        pred.row += pred_offset + 1
                        pred.col += pred_offset + 1

                        split.internal = [(pred, [])]

                        lcs_nodes[lcs_pos - 1].append(split)

                        pred = split

                    if node.outgoing == lcs_pos:
                        print(f"Split outgoing: {node}")

                        split = Node(node.row, node.col, node.len - 1)
                        split.edges = node.edges

                        split.internal = [(node, [])]
                        lcs_nodes[lcs_pos - 1].append(split)

                        pred.edges.append((node, [variant]))

                        node.row += offset
                        node.col += offset
                        node.len = 1
                        node.edges = []

                        break

                    else:

                        pred.edges.append((node, [variant]))
                        node.incoming = lcs_pos
                        pred.outgoing = lcs_pos


            for pred_idx, pred in enumerate(lcs_nodes[lcs_pos - 1]):
                pred_offset = pred.len - 1

                if node.row + offset > pred.row + pred_offset and node.col + offset > pred.col + pred_offset:

                    print(f"P{node.row + offset, node.col + offset} vs {pred.row + pred_offset, pred.col + pred_offset}")

                    variant = Variant(pred.row + pred_offset, node.row + offset - 1, observed[pred.col + pred_offset:node.col + offset - 1])
                    print(variant.to_hgvs(reference))

                    if pred.incoming == lcs_pos:
                        print(f"Split incoming {pred}")

                        split = Node(pred.row, pred.col, pred.len)

                        pred.row += pred_offset + 1
                        pred.col += pred_offset + 1

                        split.internal = [(pred, [])]

                        lcs_nodes[lcs_pos - 1][pred_idx] = split

                        pred = split

                    if node.outgoing == lcs_pos:
                        print(f"Split outgoing: {node}")

                        split = Node(node.row, node.col, node.len - 1)
                        split.edges = node.edges

                        split.internal = [(node, [])]
                        lcs_nodes[lcs_pos - 1].append(split)

                        pred.edges.append((node, [variant]))

                        node.row += offset
                        node.col += offset
                        node.len = 1
                        node.edges = []

                        break

                    else:

                        pred.edges.append((node, [variant]))
                        node.incoming = lcs_pos
                        pred.outgoing = lcs_pos

            if node.len > 1:
                node.len -= 1
                lcs_nodes[lcs_pos - 1].append(node)

            print(level,)
            for l_idx, level in enumerate(lcs_nodes):
                print(l_idx)
                for l in level:
                    print('  ', l)

    for node in lcs_nodes[0]:
        if node.edges:
            variant = Variant(0, node.row - 1, observed[:node.col - 1])
            # graph[source].append(((node.row, node.col), [variant] if variant else []))
            source.edges.append((node, [variant] if variant else []))

    return source


def traversal(reference, observed, source, atomics=False):
    def traverse(node, path):
        if not node.edges and not node.internal:
            yield path
            return

        for succ, variant in node.edges + node.internal:
            if atomics and len(variant) > 0:
                for atomic in variant[0].atomics():
                    yield from traverse(succ, path + atomic)
            else:
                yield from traverse(succ, path + variant)

    yield from traverse(source, [])


def to_dot(reference, source):
    dot = "digraph {\n"
    queue = [source]
    visited = [source]
    while queue:
        node = queue.pop(0)
        dot += f'    "{node.row}_{node.col}" [label="{node.row, node.col}"];\n'
        for succ, variant in node.edges + node.internal:
            dot += f'    "{node.row}_{node.col}" -> "{succ.row}_{succ.col}" [label="{to_hgvs(variant, reference, sequence_prefix=False)}"];\n'
            if succ not in visited:
                visited.append(succ)
                queue.append(succ)

    return dot + "}"
