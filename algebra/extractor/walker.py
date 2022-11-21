import sys
from algebra.extractor.cover import brute_cover, cover, find_pmrs, inv_array, overlapping, print_tables


def walk_sol(solutions, max_result, pos=None, path=None):

    for idx, period, cnt in solutions[pos]:

        # Adjacent entries from same pmrs are not part of minimal solution
        if path and idx == path[-1][3]:
            break

        for count in range(cnt, 1, -1):

            # Move left until we find a position with a solution
            prev_pos = pos - period * count
            while not solutions[prev_pos]:
                prev_pos -= 1
                if prev_pos < 1:
                    break

            # The actual entry
            entry = pos - period * count + 1, period, count, idx

            if prev_pos > 0:
                # There is possibly more to explore
                yield from walk_sol(solutions, max_result, prev_pos, path + [entry])
            else:
                # Otherwise, yield result if it is a maximal solution
                if sum([x[1] * x[2] for x in path + [entry]]) == max_result:
                    yield reversed(path + [entry])


def walk(inv, pmrs, max_cover, overlap, word, pos, path):

    for idx in inv[pos]:
        # Adjacent entries from same pmrs are not part of minimal solution
        if path and idx == path[-1][3]:
            continue

        begin, period, pmrs_count, _ = pmrs[idx]

        # Try all counts downwards
        for count in range(pmrs_count, 1, -1):
            start = max(pos - period * count + 1, begin)

            # Check if this candidate collides with prev
            if path and start + count * period > path[-1][0]:
                continue

            # Move left until we find a position with a solution
            prev_pos = pos - period * count
            while not inv[prev_pos]:
                prev_pos -= 1
                if prev_pos < 1:
                    break

            # The actual entry
            entry = start, period, count, idx

            if prev_pos > 0:
                # There is possibly more to explore
                yield from walk(inv, pmrs, max_cover, overlap, word, prev_pos, path + [entry])
            else:
                # Otherwise, yield result if it is a maximal solution
                if sum([x[1] * x[2] for x in path + [entry]]) == max_cover[-1]:
                    yield reversed(path + [entry])


def sol2paths(solutions, max_cover):
    for start in range(len(max_cover) - 1, 0, -1):
        if max_cover[start] != max_cover[-1]:
            # Don't look for paths without maximal result
            break
        yield from walk_sol(solutions, max_cover[-1], start, [])


def inv2paths(inv, pmrs, max_cover, overlap, word):
    for start in range(len(max_cover) - 1, 0, -1):
        if max_cover[start] != max_cover[-1]:
            # Don't look for paths without maximal result
            break
        yield from walk(inv, pmrs, max_cover, overlap, word, start, [])


def fill(subs, word):
    desc = []

    end = -1
    total = 0
    for start, period, count, _ in subs:

        # This position already covered by previous iteration
        if end >= start:
            continue

        # Add leading bases that were not part of a repeat
        if end + 1 < start:
            desc.append(word[end + 1: start])

        # Add repeat to description
        desc.append(word[start:start + period] + f"[{count}]")
        total += period * count

        # Set new end base
        end = start + period * count - 1

    # Add trailing bases that was not part of a repeat
    if end + 1 < len(word):
        desc.append(word[end + 1:])

    return desc, total


def path2hgvs(path, word):
    return ";".join(fill(path, word)[0])


def main():
    if len(sys.argv) < 2:
        print(f"usage: {sys.argv[0]} word", file=sys.stderr)
        return

    word = sys.argv[1]
    n = len(word)
    print(n, word)

    pmrs = sorted(find_pmrs(word))
    overlap = overlapping(pmrs)
    for idx, pmr in enumerate(pmrs):
        print(f"{pmr},  # {idx:2}: {word[pmr[0]:pmr[0] + pmr[1]]} : {overlap[idx]}")

    inv = inv_array(n, pmrs)
    max_cover, solutions = cover(word, pmrs, hgvs=True)
    print_tables(n, word, inv, max_cover)

    print("solutions:")
    for idx, sol in enumerate(solutions):
        print(idx, sol)
    paths = sol2paths(solutions, max_cover)
    print([path2hgvs(path, word) for path in paths])
    print("New")
    # paths = walk(inv, pmrs, max_cover, overlap, word, len(word) - 1, [])
    # print(list(paths))
    # print([path2hgvs(path, word) for path in paths])
    paths = inv2paths(inv, pmrs, max_cover, overlap, word)
    print([path2hgvs(path, word) for path in paths])


if __name__ == "__main__":
    main()
