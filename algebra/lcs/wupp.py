def edit(reference, observed):

    matrix = [[0 for _ in range(len(observed) + 1)]
              for _ in range(len(reference) + 1)]

    def func(idx):
        limit = max(diagonals[idx + offset - 1], diagonals[idx + offset + 1])
        print(f"limit: {limit}")

        if idx >= 0:
            row = diagonals[idx + offset] + 1
            col = row + idx
            print(f"row and col before while {row}, {col}")

            while (row < limit or
                   (col <= len(observed) and
                    row <= len(reference) and
                    reference[row - 1] == observed[col - 1])):
                print(row, col, reference[row - 1] == observed[col - 1])
                matrix[row][col] = delta + 2 * it
                row += 1
                col += 1

            if col <= len(observed) and row <= len(reference):
                matrix[row][col] = delta + 2 * it + 2

            return row
        else:
            col = diagonals[idx + offset] + 1
            row = col - idx

            while (col < limit or
                   (col <= len(observed) and
                    row <= len(reference) and
                    reference[row - 1] == observed[col - 1])):
                print(row, col, reference[row - 1] == observed[col - 1])
                matrix[row][col] = delta + 2 * it
                row += 1
                col += 1

            if col <= len(observed) and row <= len(reference):
                matrix[row][col] = delta + 2 * it + 2

            return col

    diagonals = [0] * (len(reference) + len(observed) + 3)
    offset = len(reference) + 1
    delta = len(observed) - len(reference)

    it = 0

    while diagonals[delta + offset] < len(observed):
        print(f"it: {it}")

        for diag_idx in range(-it, delta):
            print(f"lower: {diag_idx}")
            diagonals[diag_idx + offset] = func(diag_idx)

        for diag_idx in range(delta + it, delta, -1):
            print(f"upper: {diag_idx}")
            diagonals[diag_idx + offset] = func(diag_idx)

        print(f"delta: {delta}")
        diagonals[delta + offset] = func(delta)

        it += 1

    from pprint import pprint
    pprint(matrix)
    print(diagonals)