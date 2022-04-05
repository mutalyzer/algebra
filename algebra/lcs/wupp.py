def edit(reference, observed):
    def extend(idx):
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

        print(f"    start at {row, col}")
        for _ in range(start, end):
            print(f"        implicit {row, col}")
            matrix[row + 1][col + 1] = abs(delta) + 2 * it
            row += 1
            col += 1

        print(f"    continue with {row, col}")
        steps = end + 1
        while row < len(reference) and col < len(observed) and reference[row] == observed[col]:
            print(f"        match {row, col}")
            matrix[row + 1][col + 1] = abs(delta) + 2 * it
            row += 1
            col += 1
            steps += 1

        print(f"    last checked {row, col}")
        matrix[row + 1][col + 1] = abs(delta) + 2 * it + 2
        return steps

    matrix = [[None for _ in range(len(observed) + 2)] for _ in range(len(reference) + 2)]

    delta = len(observed) - len(reference)
    offset = len(reference) + 1
    diagonals = [0] * (len(reference) + len(observed) + 3)
    it = 0

    while diagonals[offset + delta] <= max(len(reference), len(observed)) - abs(delta):

        if delta >= 0:
            lower = range(-it, delta)
            upper = range(delta + it, delta, -1)
        else:
            lower = range(delta - it, delta)
            upper = range(it, delta, -1)

        for idx in lower:
            print(f"extend diagonal {idx} (lower)")
            diagonals[offset + idx] = extend(idx)
        for idx in upper:
            print(f"extend diagonal {idx} (upper)")
            diagonals[offset + idx] = extend(idx)
        print(f"extend diagonal {delta} (delta)")
        diagonals[offset + delta] = extend(delta)

        for row in matrix:
            print(row)
        print(diagonals)

        it += 1

    return abs(delta) + 2 * (it - 1), [r[:-1] for r in matrix[:-1]]
