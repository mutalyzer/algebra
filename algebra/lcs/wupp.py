def edit(reference, observed):
    def extend(idx, start, end, value):
        if idx >= 0:
            row = start
            col = row + idx
        else:
            col = start
            row = col - idx

        steps = 1
        print(f"    start at {row, col}")
        for _ in range(start, end):
            print(f"        implicit {row, col}")
            matrix[row][col] = value
            row += 1
            col += 1
            steps += 1

        print(f"    continue with {row, col}")
        while row <= len(reference) and col <= len(observed) and reference[row - 1] == observed[col - 1]:
            print(f"        match {row, col}")
            matrix[row][col] = value
            row += 1
            col += 1
            steps += 1

        print(f"    last checked {row, col}")
        #matrix[row][col] = value + 2
        return start + steps

    matrix = [[None for _ in range(len(observed) + 1)] for _ in range(len(reference) + 1)]

    delta = len(observed) - len(reference)
    offset = len(reference) + 1
    diagonals = [1] * (len(reference) + len(observed) + 3)
    it = 0

    while diagonals[offset + delta] <= max(len(reference), len(observed)):
        for idx in range(-it, delta):
            print(f"extend diagonal {idx} (lower)")
            diagonals[offset + idx] = extend(idx, diagonals[offset + idx], max(diagonals[offset + idx - 1] + 1, diagonals[offset + idx + 1]) - 1, abs(delta) + 2 * it)
        for idx in range(delta + it, delta, -1):
            print(f"extend diagonal {idx} (upper)")
            diagonals[offset + idx] = extend(idx, diagonals[offset + idx], max(diagonals[offset + idx - 1], diagonals[offset + idx + 1] + 1) - 1, abs(delta) + 2 * it)
        print(f"extend diagonal {delta} (delta)")
        diagonals[offset + delta] = extend(delta, diagonals[offset + delta], max(diagonals[offset + delta - 1], diagonals[offset + delta + 1]) - delta, abs(delta) + 2 * it)

        for row in matrix:
            print(row)
        print(diagonals)

        it += 1

    return abs(delta) + 2 * (it - 1), matrix
