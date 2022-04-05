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

        steps = 1
        print(f"    start at {row, col}")
        for _ in range(start, end):
            print(f"        implicit {row, col}")
            matrix[row][col] = abs(delta) + 2 * it
            row += 1
            col += 1
            steps += 1

        print(f"    continue with {row, col}")
        while row <= len(reference) and col <= len(observed) and reference[row - 1] == observed[col - 1]:
            print(f"        match {row, col}")
            matrix[row][col] = abs(delta) + 2 * it
            row += 1
            col += 1
            steps += 1

        print(f"    last checked {row, col}")
        matrix[row][col] = abs(delta) + 2 * it + 2
        return start + steps

    matrix = [[None for _ in range(len(observed) + 2)] for _ in range(len(reference) + 2)]

    delta = len(observed) - len(reference)
    offset = len(reference) + 1
    diagonals = [1] * (len(reference) + len(observed) + 3)
    it = 0

    while diagonals[offset + delta] <= len(observed) + 1 - delta:
        for idx in range(-it, delta):
            print(f"extend diagonal {idx} (lower)")
            diagonals[offset + idx] = extend(idx)
        for idx in range(delta + it, delta, -1):
            print(f"extend diagonal {idx} (upper)")
            diagonals[offset + idx] = extend(idx)
        print(f"extend diagonal {delta} (delta)")
        diagonals[offset + delta] = extend(delta)

        for row in matrix:
            print(row)
        print(diagonals)

        it += 1

    return abs(delta) + 2 * (it - 1), [r[:-1] for r in matrix[:-1]]
