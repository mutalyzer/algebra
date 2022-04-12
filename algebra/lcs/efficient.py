from ..variants.variant import Variant


class Node:
    def __init__(self, row, col):
        self.row = row
        self.col = col
        self.child = []

    def __str__(self):
        return f"{self.row, self.col}"


def heur_func(reference, observed, row, col):
    """Heuristic function for A*.
    Args:
        reference (string): reference sequence
        observed (string): observed sequence
        row (integer): row of the element
        col (integer): column of the element
    Returns: the minimal edit distance for the element at (row, col) to
             the end of both (reference and observed) strings.
    """
    return abs((len(reference) - row) - (len(observed) - col))


def edit(reference, observed):
    """Calculates all longest common subsequences (LCS).
    Args:
        reference (string): reference sequence
        observed (string): observed sequence
    """
    f = heur_func(reference, observed, 0, 0)

    graph = [[] for _ in range(min(len(reference), len(observed)))]
    lcs_len = 0

    rows = [1] + [-1] * len(reference)
    cols = [1] + [-1] * len(observed)

    row_start = 0
    col_start = 0
    row_end = 1
    col_end = 1

    offset = 0

    while True:
        for row in range(row_start, row_end):
            col = rows[row] + offset

            if col > len(observed):
                if row == row_start:
                    row_start = row + 1
                continue

            if row > 0 and reference[row - 1] == observed[col - 1]:
                if rows[row - 1] + offset > col - 1:
                    distance = f - heur_func(reference, observed, row, col)
                    lcs_pos = ((row + col) - distance) // 2 - 1
                    graph[lcs_pos].append(Node(row, col))
                    lcs_len = max(lcs_len, lcs_pos + 1)
                elif row == row_start:
                    distance = (f + 2) - heur_func(reference, observed, row, col)
                    lcs_pos = ((row + col) - distance) // 2 - 1
                    graph[lcs_pos].append(Node(row, col))
                    lcs_len = max(lcs_len, lcs_pos + 1)

            if (
                heur_func(reference, observed, row, col)
                < heur_func(reference, observed, row, col - 1)
                or (
                    row > 0
                    and rows[row - 1] + offset > col
                    and heur_func(reference, observed, row, col)
                    < heur_func(reference, observed, row - 1, col)
                )
                or (
                    row > 0
                    and rows[row - 1] + offset > col - 1
                    and reference[row - 1] == observed[col - 1]
                )
            ):
                rows[row] += 1
                cols[col] = max(cols[col] + offset, row + 1) - offset
                col_end = max(col_end, col + 1)

            elif row == row_start:
                row_start = row + 1

        for col in range(col_start, col_end):
            row = cols[col] + offset

            if row > len(reference):
                if col == col_start:
                    col_start = col + 1
                continue

            if col > 0 and reference[row - 1] == observed[col - 1]:
                if cols[col - 1] + offset > row - 1:
                    distance = f - heur_func(reference, observed, row, col)
                    lcs_pos = ((row + col) - distance) // 2 - 1
                    graph[lcs_pos].append(Node(row, col))
                    lcs_len = max(lcs_len, lcs_pos + 1)
                elif col == col_start:
                    distance = (f + 2) - heur_func(reference, observed, row, col)
                    lcs_pos = ((row + col) - distance) // 2 - 1
                    graph[lcs_pos].append(Node(row, col))
                    lcs_len = max(lcs_len, lcs_pos + 1)

            if (
                heur_func(reference, observed, row, col)
                < heur_func(reference, observed, row - 1, col)
                or (
                    col > 0
                    and cols[col - 1] + offset > row
                    and heur_func(reference, observed, row, col)
                    < heur_func(reference, observed, row, col - 1)
                )
                or (
                    col > 0
                    and cols[col - 1] + offset > row - 1
                    and reference[row - 1] == observed[col - 1]
                )
            ):
                cols[col] += 1
                rows[row] = max(rows[row] + offset, col + 1) - offset
                row_end = max(row_end, row + 1)

            elif col == col_start:
                col_start = col + 1

        if row_start > len(reference) and col_start > len(observed):
            break

        if row_start == row_end and col_start == col_end:
            if reference[row_end - 1] == observed[col_end - 1]:
                lcs_pos = (
                    (row_end + col_end)
                    - (f - heur_func(reference, observed, row_end, col_end))
                ) // 2 - 1
                graph[lcs_pos].append(Node(row_end, col_end))
                lcs_len = max(lcs_len, lcs_pos + 1)

                row_end += 1
                col_end += 1
                rows[row_start] = col_end - offset
                cols[col_start] = row_end - offset
                continue

            f += 2
            offset += 1
            rows[row_end] = col_end - offset
            cols[col_end] = row_end - offset

            row_start = 0
            col_start = 0
            row_end += 1
            col_end += 1

    return f, graph[:lcs_len]


def manhattan(observed, row_start, col_start, row_end, col_end):
    ops = set()

    for row in range(row_start + 1, row_end):
        ops.add((row, "del"))

    for row in range(row_start, row_end):
        for col in range(col_start, col_end - 1):
            ops.add((row, "ins", observed[col]))

    return ops


def build(graph, reference, observed):
    ops = set()

    if graph == []:
        return manhattan(observed, 0, 0, len(reference) + 1, len(observed) + 1), graph

    # all nodes on the last level are reachable and they indicate the
    # end of the LCS
    for node in graph[len(graph) - 1]:
        ops.update(
            manhattan(
                observed, node.row, node.col, len(reference) + 1, len(observed) + 1
            )
        )
        node.child = [None]

    # construct the LCSs from the last match to the first match
    for level in range(len(graph), 1, -1):
        # unPythonic allows for removal of nodes while iterating
        idx = 0
        while idx < len(graph[level - 1]):
            node = graph[level - 1][idx]
            if len(node.child) > 0:
                for prev in graph[level - 2]:
                    if prev.row < node.row and prev.col < node.col:
                        ops.update(
                            manhattan(observed, prev.row, prev.col, node.row, node.col)
                        )
                        prev.child.append(node)
                idx += 1
            else:
                # remove unused node
                graph[level - 1].pop(idx)

    # remove unused nodes on root level
    graph[0] = [node for node in graph[0] if len(node.child) != 0]

    for node in graph[0]:
        ops.update(manhattan(observed, 0, 0, node.row, node.col))

    return ops, graph


def traversal(reference, observed, graph, atomics=False):
    def traverse(node, path):
        for child in node.child:
            if child is None:
                variant = Variant(node.row, len(reference), observed[node.col:])
                if variant:
                    if atomics:
                        for atomic in variant.atomics():
                            yield path + atomic
                    else:
                        yield path + [variant]
                else:
                    yield path
                return

            variant = Variant(node.row, child.row - 1, observed[node.col:child.col - 1])
            if variant:
                if atomics:
                    for atomic in variant.atomics():
                        yield from traverse(child, path + atomic)
                else:
                    yield from traverse(child, path + [variant])
            else:
                yield from traverse(child, path)

    if graph == []:
        variant = Variant(0, len(reference), observed)
        if atomics:
            for atomic in variant.atomics():
                yield atomic
        else:
            yield [variant]
        return

    for node in graph[0]:
        variant = Variant(0, node.row - 1, observed[:node.col - 1])
        if atomics and variant:
            for atomic in variant.atomics():
                yield from traverse(node, atomic)
        elif variant:
            yield from traverse(node, [variant])
        else:
            yield from traverse(node, [])
