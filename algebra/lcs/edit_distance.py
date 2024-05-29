"""Calculate the simple edit distance.

Calculates the simple edit distance by calculating a Longest Common
Subsequence (LCS) alignment [1]_. The implementation is adapted from [2]_.
This method is useful when only the simple edit distance is needed (not
an or all alignments).

See Also
--------
algebra.lcs.lcs_graph : Calculate all LCS alignments.

References
----------
[1] S. Wu, U. Manber, G. Myers and W. Miller. "An O(NP) Sequence
Comparison Algorithm". In: Information Processing Letters, 35.6 (1990),
pp. 317-323.
[2] https://github.com/cubicdaiya/onp/blob/master/python/onp.py.
"""


def edit_distance(a, b):
    """Calculate the simple edit distance between two strings."""
    def snake(k, lower, upper):
        col = max(lower, upper)
        row = col - k
        while row < m and col < n and a[row] == b[col]:
            row += 1
            col += 1
        return col

    m = len(a)
    n = len(b)
    if m >= n:
        a, b = b, a
        m, n = n, m

    offset = m + 1
    delta = n - m
    fp = [-1] * (m + n + 3)
    p = -1

    while True:
        p += 1
        for k in range(-p, delta):
            fp[k + offset] = snake(k, fp[k - 1 + offset] + 1, fp[k + 1 + offset])
        for k in range(delta + p, delta, -1):
            fp[k + offset] = snake(k, fp[k - 1 + offset] + 1, fp[k + 1 + offset])
        fp[delta + offset] = snake(delta, fp[delta - 1 + offset] + 1, fp[delta + 1 + offset])

        if fp[delta + offset] >= n:
            break

    return delta + 2 * p
