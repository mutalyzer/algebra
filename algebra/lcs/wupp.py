from ..variants.variant import Variant, to_hgvs


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
    graph = {sink: []}

    source = (0, 0)
    if lcs_nodes == [[]]:
        graph[source] = [(sink, [Variant(0, len(reference), observed)])]
        return graph

    for node in lcs_nodes[-1]:
        variant = Variant(node["row"] + node["len"] - 1, len(reference), observed[node["col"] + node["len"] - 1:])
        graph[(node["row"], node["col"])] = [(sink, [variant] if variant else [])]

    print(to_dot(reference, graph))

    graph[source] = []
    for idx, nodes in enumerate(lcs_nodes[:0:-1]):
        level = len(lcs_nodes) - idx - 1
        print(f"Entering level: {level}")

        for node in list(nodes):
            if (node["row"], node["col"]) not in graph:
                print(f"REMOVE {node}")
                nodes.remove(node)

        for node in nodes:
            for target in lcs_nodes[level - 1]:
                offset = node["len"] - (node["lcs_pos"] - level) - 1
                target_offset = target["len"] - (target["lcs_pos"] - level + 1) - 1

                if node["row"] + offset <= target["row"] + target_offset + 1 and node["col"] + offset <= target["col"] + target_offset + 1:
                    continue

                print(f"{node['row'] + offset, node['col'] + offset} -> {target['row'] + target_offset, target['col'] + target_offset}")

                dest = (node['row'] + offset, node['col'] + offset)
                if target_offset < target["len"] - 1:
                    print("SPLIT target")
                    split = (target["row"] + target_offset + 1, target["col"] + target_offset + 1)
                    print(split)
                    target["lcs_pos"] -= target_offset
                    target["len"] -= target_offset
                    print(target)
                    graph[split] = graph[(target["row"], target["col"])]
                    graph[(target["row"], target["col"])] = [(split, [])]

                if offset > 0:
                    print("SPLIT node")
                    split = (node["row"] + offset, node["col"] + offset)
                    print(split)
                    node["lcs_pos"] -= offset - 1
                    node["len"] -= offset - 1
                    print(node)
                    graph[split] = graph[(node["row"], node["col"])]
                    graph[(node["row"], node["col"])] = [(split, [])]
                    offset = 0
                    dest = split

                print("EDGE")
                variant = Variant(target["row"] + target_offset, dest[0] + offset - 1, observed[target["col"] + target_offset:dest[1] + offset - 1])
                if (target["row"], target["col"]) not in graph:
                    graph[(target["row"], target["col"])] = []
                graph[(target["row"], target["col"])].append((dest, [variant]))

                print(to_dot(reference, graph))

    for node in lcs_nodes[0]:
        if (node["row"], node["col"]) in graph:
            variant = Variant(0, node["row"] - 1, observed[:node["col"] - 1])
            graph[source].append(((node["row"], node["col"]), [variant] if variant else []))

    print(to_dot(reference, graph))

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


def to_dot(reference, graph):
    dot = "digraph {\n"
    for node, edges in graph.items():
        dot += f'    "{node[0]}_{node[1]}" [label="{node}"];\n'
        for child, edge in edges:
            dot += f'    "{node[0]}_{node[1]}" -> "{child[0]}_{child[1]}" [label="{to_hgvs(edge, reference, sequence_prefix=False)}"];\n'
    return dot + "}"
