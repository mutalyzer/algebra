def edit(reference, observed):
    """Longest Common Subsequence (LCS) edit distance
    The traditional dynamic programming approach

        Parameters:
            reference (str):
            observed (str):

        Returns:
            The edit distance matrix
    """
    d = [[None for _ in range(len(observed) + 1)] for _ in range(len(reference) + 1)]

    for i in range(len(reference) + 1):
        d[i][0] = i

    for i in range(len(observed) + 1):
        d[0][i] = i

    for i in range(len(reference)):
        for j in range(len(observed)):
            if reference[i] == observed[j]:
                d[i + 1][j + 1] = min(d[i][j + 1] + 1, d[i + 1][j] + 1, d[i][j])
            else:
                d[i + 1][j + 1] = min(d[i][j + 1] + 1, d[i + 1][j] + 1)

    return d


def traversal(reference, observed, d):
    pass
