from ..variants.variant import Variant


def edit(reference, observed):
    """Longest Common Subsequence (LCS) edit distance
    """
    def heuristic(row, col):
        """Heuristic function for A*
        """
        return abs((len(reference) - row) - (len(observed) - col))

    def lcs_match(lcs_nodes, lcs_len, row, col, f):
        """
        """
        distance = f - heuristic(row, col)
        lcs_pos = ((row + col) - distance) // 2 - 1
        length = lcs_nodes[lcs_pos - 1].pop((row - 1, col - 1), {"length": 0})["length"]
        lcs_nodes[lcs_pos][(row, col)] = {"length": length + 1}
        return max(lcs_len, lcs_pos + 1)

    f = heuristic(0, 0)

    lcs_nodes = [{} for _ in range(min(len(reference), len(observed)))]
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
                    lcs_len = lcs_match(lcs_nodes, lcs_len, row, col, f)
                elif row == row_start:
                    lcs_len = lcs_match(lcs_nodes, lcs_len, row, col, f + 2)

            if (heuristic(row, col) < heuristic(row, col - 1) or
                    (row > 0 and rows[row - 1] + offset > col and heuristic(row, col) < heuristic(row - 1, col)) or
                    (row > 0 and rows[row - 1] + offset > col - 1 and reference[row - 1] == observed[col - 1])):
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
                    lcs_len = lcs_match(lcs_nodes, lcs_len, row, col, f)
                elif col == col_start:
                    lcs_len = lcs_match(lcs_nodes, lcs_len, row, col, f + 2)

            if (heuristic(row, col) < heuristic(row - 1, col) or
                    (col > 0 and cols[col - 1] + offset > row and heuristic(row, col) < heuristic(row, col - 1)) or
                    (col > 0 and cols[col - 1] + offset > row - 1 and reference[row - 1] == observed[col - 1])):
                cols[col] += 1
                rows[row] = max(rows[row] + offset, col + 1) - offset
                row_end = max(row_end, row + 1)

            elif col == col_start:
                col_start = col + 1

        if row_start > len(reference) and col_start > len(observed):
            break

        if row_start == row_end and col_start == col_end:
            if reference[row_end - 1] == observed[col_end - 1]:
                lcs_len = lcs_match(lcs_nodes, lcs_len, row_end, col_end, f)

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

    return f, lcs_nodes[:lcs_len]


def lcs_graph(reference, observed, lcs_nodes):
    if lcs_nodes == []:
        return {}

    graph = {(len(reference) + 1, len(observed) + 1): list(lcs_nodes[-1].keys())}
    for level in range(len(lcs_nodes) - 1, 0, -1):
        for to_node, to_length in lcs_nodes[level].items():
            for length in range(to_length):
                for from_node, from_length in lcs_nodes[level - length - 1].items():
                    if from_node[0] < to_node[0] - length and from_node[1] < to_node[1] - length:
                        print(f"edge {from_node}:{from_length} -> {to_node}:{to_length}")
                        from_node = True

    print(lcs_nodes)
    return graph
