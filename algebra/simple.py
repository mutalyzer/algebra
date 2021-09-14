import sys


def subsequence(list1, list2):
    """The length of the longest prefix of list1 in list2."""

    matched = []
    start = 0
    for item in list1:
        if start >= len(list2):
            return matched
        try:
            start = list2.index(item, start) + 1
        except ValueError:
            continue
        matched.append(item)
    return matched


def sol_compare(observed1, solutions1, observed2, solutions2):
    if observed1 == observed2:
        return "equivalent"
    # empty paths are a special case, because the empty path is always a
    # subpath of any path, but not something we are interessed in
    if solutions1 == [[]] or solutions2 == [[]]:
        return "disjoint"

    overlap = []
    for path1 in solutions1:
        for path2 in solutions2:
            matched = subsequence(path1, path2)
            if len(path1) == len(matched):
                return "is_contained"

            if len(matched) > len(overlap):
                overlap = matched

            matched = subsequence(path2, path1)
            if len(path2) == len(matched):
                return "contains"

            if len(matched) > len(overlap):
                overlap = matched

    if len(overlap) > 0:
        return "overlap"
    return "disjoint"


def edit(ref, obs):
    """Longest common subsequence edit distance.
    The traditional dynamic programming solution.
    """
    d = [[None for i in range(len(ref) + 1)] for j in range(len(obs) + 1)]

    for i in range(len(ref) + 1):
        d[0][i] = i

    for i in range(len(obs) + 1):
        d[i][0] = i

    for i in range(len(obs)):
        for j in range(len(ref)):
            if ref[j] == obs[i]:
                d[i + 1][j + 1] = min(d[i][j + 1] + 1, d[i + 1][j] + 1, d[i][j])
            else:
                d[i + 1][j + 1] = min(d[i][j + 1] + 1, d[i + 1][j] + 1)
    return d


def traverse(ref, obs, d, i, j, path, sols):
    """Recursively enumerate all minimal edit scripts."""
    if i == 0 and j == 0:
        sols.append(path.copy())
        return 1

    count = 0
    if i > 0 and j > 0 and d[i][j] == d[i - 1][j - 1] and ref[j - 1] == obs[i - 1]:
        count += traverse(ref, obs, d, i - 1, j - 1, path, sols)

    if i > 0 and d[i][j] == d[i - 1][j] + 1:
        path.append((j, "ins", obs[i - 1]))
        count += traverse(ref, obs, d, i - 1, j, path, sols)
        path.pop()

    if j > 0 and d[i][j] == d[i][j - 1] + 1:
        path.append((j, "del"))
        count += traverse(ref, obs, d, i, j - 1, path, sols)
        path.pop()

    return count


def compare(reference, lhs, rhs, debug=False):
    lhs_matrix = edit(reference, lhs)
    lhs_paths = []
    traverse(reference, lhs, lhs_matrix, len(lhs), len(reference), [], lhs_paths)
    a = []
    for x in lhs_paths:
        for y in x:
            a.append(y)
    lhs_elements = set(a)
    if debug:
        print(f"lhs_elements ({len(lhs_elements)}): {lhs_elements}")
        print(f"lhs_paths ({len(lhs_paths)}): {lhs_paths}")

    rhs_matrix = edit(reference, rhs)
    rhs_paths = []
    traverse(reference, rhs, rhs_matrix, len(rhs), len(reference), [], rhs_paths)
    a = []
    for x in rhs_paths:
        for y in x:
            a.append(y)
    rhs_elements = set(a)
    if debug:
        print(f"rhs_elements ({len(rhs_elements)}): {rhs_elements}")
        print(f"rhs_paths ({len(rhs_paths)}): {rhs_paths}")

    return sol_compare(lhs, lhs_paths, rhs, rhs_paths), lhs_elements, rhs_elements


if __name__ == "__main__":
    print(compare(sys.argv[1], sys.argv[2], sys.argv[3]))
