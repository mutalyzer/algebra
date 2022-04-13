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
                node = {'row': match_row + 1, 'col': match_col + 1, 'len': row - match_row, 'lcs_pos': lcs_pos}
                for i in range(row - match_row):
                    lcs_nodes[lcs_pos - i].append(node)
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
            node = {'row': match_row + 1, 'col': match_col + 1, 'len': row - match_row, 'lcs_pos': lcs_pos}
            for i in range(row - match_row):
                lcs_nodes[lcs_pos - i].append(node)

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
        variant = Variant(node["row"] + node["len"] - 1, len(reference), observed[node["col"] + node["len"] - 1:])
        graph[(node["row"], node["col"])] = [(sink, [variant] if variant else [])]
        print(variant.to_hgvs(reference))

    for idx, nodes in enumerate(lcs_nodes[::-1]):
        level = len(lcs_nodes) - idx - 1
        print(f"Entering level: {level}")

        for node in nodes[:]:
            if (node["row"], node["col"]) not in graph:
                print(f"REMOVE {node}")
                nodes.remove(node)
                continue

            print(f"Node: {node}")

            if level == 0:
                variant = Variant(0, node["row"] - 1, observed[:node["col"] - 1])
                graph[(0, 0)].append(((node["row"], node["col"]), [variant] if variant else []))
                print(variant.to_hgvs(reference))
            else:
                offset = node["len"] - (node["lcs_pos"] - level) - 1
                print(f"    node {node['row'] + offset, node['col'] + offset} @ offset {offset}")

                for target in lcs_nodes[level - 1]:
                    # Skip self
                    if node is target:
                        continue

                    target_offset = target["len"] - (target["lcs_pos"] - level + 1) - 1
                    print(f"        target {target['row'] + target_offset, target['col'] + target_offset} @ offset {target_offset}")

                    node_row = node["row"] + offset
                    target_row = target["row"] + target_offset

                    node_col = node["col"] + offset
                    target_col = target["col"] + target_offset

                    if node_row > target_row and node_col > target_col:
                        target_coor = target["row"], target["col"]
                        # target_coor = target_row, target_col
                        if target_coor not in graph:
                            graph[target_coor] = []
                        variant = Variant(target_row, node_row - 1, observed[target_col:node_col - 1])
                        graph[target_coor].append(((node["row"], node["col"]), [variant]))
                        # graph[target_coor].append(((node_row, node_col), [variant]))
                        print(target_coor)
                        print(variant.to_hgvs(reference))

    return graph


def traversal(reference, observed, graph, atomics=False):
    def traverse(node, path):
        if node == (len(reference) + 1, len(observed) + 1):
            # print("final", node)
            yield path
            return

        for child, variant in graph[node]:
            # print("child", child)
            if atomics and len(variant) > 0:
                for atomic in variant[0].atomics():
                    yield from traverse(child, path + atomic)
            else:
                yield from traverse(child, path + variant)

    yield from traverse((0, 0), [])
