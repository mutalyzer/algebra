def edit(a, b):
    """Longest Common Subsequence (LCS) edit distance
    Based on 'An O(NP) Sequence Comparison Algorithm' by Sun Wu, Udi Manber and Gene Myers
    """
    def snake(m, n, k, pp, ppp):
        y = max(pp, ppp)
        x = y - k
        while x < m and y < n and a[x] == b[y]:
            x += 1
            y += 1
        return y

    m = len(a)
    n = len(b)
    if m >= n:
        a, b = b, a
        m, n = n, m

    offset = m + 1
    delta = n - m
    size = m + n + 3
    fp = [-1] * size
    p = -1

    while True:
        p += 1
        for k in range(-p, delta):
            fp[k + offset] = snake(m, n, k, fp[k - 1 + offset] + 1, fp[k + 1 + offset])
        for k in range(delta + p, delta, -1):
            fp[k + offset] = snake(m, n, k, fp[k - 1 + offset] + 1, fp[k + 1 + offset])
        fp[delta + offset] = snake(m, n, delta, fp[delta - 1 + offset] + 1, fp[delta + 1 + offset])

        if fp[delta + offset] >= n:
            break

    return delta + 2 * p
