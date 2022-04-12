from ..variants.variant import Variant


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
                lcs_nodes[lcs_pos].append({'row': match_row + 1, 'col': match_col + 1, 'len': row - match_row})
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
            lcs_nodes[lcs_pos].append({'row': match_row + 1, 'col': match_col + 1, 'len': row - match_row})

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

    sink = (len(reference) + 1, len(observed) + 1)
    graph = {sink: [], (0, 0): []}

    if lcs_nodes == [[]]:
        graph[(0, 0)].append((sink, [Variant(0, len(reference), observed)]))

    for node in lcs_nodes[-1]:
        variant = Variant(node['row'] + node['len'] - 1, len(reference), observed[node['col'] + node['len'] - 1:])
        graph[(node["row"], node["col"])] = [(sink, [variant] if variant else [])]

    for idx, nodes in enumerate(lcs_nodes[::-1]):
        level = len(lcs_nodes) - idx - 1
        print(f"Entering level: {level}")
        print()
        for node in nodes[:]:
            if (node["row"], node["col"]) not in graph:
                nodes.remove(node)
                continue

            print(f'Node: {node}')
            for offset in range(node['len']):
                child_row = node['row'] + node['len'] - 1 - offset
                child_col = node['col'] + node['len'] - 1 - offset
                print(f'offset: {offset} {child_row, child_col}')
                max_tgt_lvl = level
                min_tgt_lvl = max(level - offset - 1, 0)
                print(f"    min/max target level: {min_tgt_lvl}/{max_tgt_lvl}")

                if (level - offset - 1) < 0:
                    variant = Variant(0, child_row - 1, observed[:child_col - 1])
                    graph[(0, 0)].append(((node["row"], node["col"]), [variant] if variant else []))

                for tgt_level in range(max_tgt_lvl, min_tgt_lvl - 1, -1):
                    print(f"    Target level: {tgt_level}")
                    for tgt_node in lcs_nodes[tgt_level]:
                        if node == tgt_node:
                            # Skip self
                            continue
                        print('         Target node:', tgt_node)

                        tgt_offset = tgt_level - min_tgt_lvl
                        if tgt_offset >= tgt_node["len"]:
                            continue

                        tgt_row = tgt_node['row'] + tgt_node['len'] - 1 - tgt_offset
                        tgt_col = tgt_node['col'] + tgt_node['len'] - 1 - tgt_offset
                        print(f'            Target offset: {tgt_offset} level: {tgt_level - tgt_offset} {tgt_row, tgt_col}')

                        if child_row > tgt_row and child_col > tgt_col:
                            tgt_coor = tgt_node["row"], tgt_node["col"]
                            if tgt_coor not in graph:
                                graph[tgt_coor] = []
                            variant = Variant(tgt_row, child_row - 1, observed[tgt_col:child_col - 1])
                            graph[tgt_coor].append(((node["row"], node["col"]), [variant]))
                            print(variant)
            print()

    return graph


def traversal(reference, observed, graph, atomics=False):
    def traverse(node, path):
        if node == (len(reference) + 1, len(observed) + 1):
            yield path
            return

        for child, variant in graph[node]:
            if atomics and len(variant) > 0:
                for atomic in variant[0].atomics():
                    yield from traverse(child, path + atomic)
            else:
                yield from traverse(child, path + variant)

    yield from traverse((0, 0), [])
