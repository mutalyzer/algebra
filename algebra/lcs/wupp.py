from algebra.variants.variant import Variant


def edit(reference, observed):
    def lcs_idx(row, col):
        return ((row + col) - (abs(delta) + 2 * it - abs((len(reference) - row) - (len(observed) - col)))) // 2 - 1

    def extend(idx):
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
            diagonals[offset + idx] = extend(idx)
        for idx in range(upper + it, delta, - 1):
            # print(f"extend diagonal {idx} (upper)")
            diagonals[offset + idx] = extend(idx)
        # print(f"extend diagonal {delta} (delta)")
        diagonals[offset + delta] = extend(delta)

        # for row in matrix:
        #     print(row)
        # print(diagonals)

        it += 1

    return abs(delta) + 2 * (it - 1), [r[:-1] for r in matrix[:-1]], lcs_nodes[:max_lcs_pos + 1]


def lcs_graph(reference, observed, lcs_nodes):

    def connect(child, potential):
        print(f"pair: {child}, {potential}")

        for child_offset in range(child['len']):
            for pot_offset in range(potential['len'] - 1, -1, -1):
                child_row = child['row'] + child_offset
                child_col = child['col'] + child_offset
                pot_row = potential['row'] + pot_offset
                pot_col = potential['col'] + pot_offset
                print(f"child offset: {child_offset} at {child_row, child_col}, potential offset: {pot_offset} at {pot_row, pot_col}")

                if child_row > pot_row and child_col > pot_col:
                    print(f"child dominates")

                    var = Variant(child_row, len(reference), observed[child_col:child_col + child_offset]).to_hgvs(reference)

                    if 'children' not in potential:
                        potential['children'] = []
                    else:
                        print('parent already has children')
                    potential['children'].append((child['row'], child['col']))
                    break
                else:
                    print("child does not dominate")
            else:
                # only executed if the inner loop did not break
                continue
            # only executed if the inner loop DID break
            # if child_offset == 0:
            #     print("child offset is 0, break")
            #     break


    for node in lcs_nodes[-1]:
        print(node)
        variant = Variant(node['row'] + node['len'] - 1, len(reference), observed[node['col'] + node['len'] - 1:]).to_hgvs(reference)
        print(variant)
        node['children'] = [('Sink', variant)]

    print()

    for idx, nodes in enumerate(lcs_nodes[::-1]):
        level = len(lcs_nodes) - idx - 1
        print(f"Entering level: {level}")
        for node in nodes:
            # if 'children' not in node:
            # TODO: remove?!
            #     print(f"Node {node} has no children")
            #     continue

            print(f'Node: {node}')
            for offset in range(node['len'] + 1):
                haystack_level = level - offset
                print(f'Target level {haystack_level}')
                for haystack in lcs_nodes[haystack_level]:
                    if node == haystack:
                        # Skip self
                        continue
                    connect(node, haystack)
            print()

    for level in lcs_nodes:
        print(level)
