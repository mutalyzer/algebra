import sys
from algebra.extractor.cover import cover, find_pmrs, inv_array, overlapping, print_tables, cover_length


def walk(inv, pmrs, overlap, word, pos, path):
    # print(f"pos: {pos}")

    if pos < 1:
        yield list(reversed(path))
        return

    yield from walk(inv, pmrs, overlap, word, pos - 1, path)

    # print(f"path: {path}")
    for idx in inv[pos]:

        # Adjacent entries from same pmrs are not part of minimal solution
        if path and idx == path[-1][3]:
            # print("skip because of same pmrs")
            continue

        begin, period, _, _ = pmrs[idx]

        # Don't let candidate collide with previous entry
        max_count = (pos - begin + 1) // period

        # Try all counts downwards
        for count in range(max_count, 1, -1):
            # print(f"pos: {pos} count: {count}")

            # Move left until we find a position with a solution
            prev_pos = pos - period * count
            while not inv[prev_pos]:
                prev_pos -= 1
                # print(f"shifting left: {prev_pos}")
                if prev_pos < 1:
                    break

            # The actual entry
            s = word[pos - period * count + 1: pos + 1]
            entry = pos - period * count + 1, period, count, idx, s
            # print(f"entry: {entry} ({s})")

            yield from walk(inv, pmrs, overlap, word, prev_pos, path + [entry])


def inv2paths(inv, pmrs, max_cover, overlap, word):
    for start in range(len(max_cover) - 1, 0, -1):
        if max_cover[start] != max_cover[-1]:
            # Don't look for paths without maximal result
            break
        yield from walk(inv, pmrs, overlap, word, start, [])


def fill(subs, word):
    desc = []

    end = -1
    total = 0
    for start, period, count, *_ in subs:

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


def path2gaps(path, word):
    gaps = 0
    end = -1
    for start, period, count, _ in path:

        # This position already covered by previous iteration
        if end >= start:
            continue

        # Add leading bases that were not part of a repeat
        if end + 1 < start:
            gaps += 1

        # Set new end base
        end = start + period * count - 1

    # Add trailing bases that was not part of a repeat
    if end + 1 < len(word):
        gaps += 1

    return gaps


def path2hgvs(path, word):
    return ";".join(fill(path, word)[0])


def metrics(paths, word, pmrs):
    max_repeats = 0
    min_repeats = len(pmrs)

    min_count = len(word)
    max_count = 0

    min_period = len(word)
    max_period = 0

    for path in paths:
        total = 0
        p = list(path)
        print(path2hgvs(p, word))
        # print("Number of repeat units:", len(p))

        max_repeats = max(max_repeats, len(p))
        min_repeats = min(min_repeats, len(p))

        min_entry_count = len(word)
        max_entry_count = 0

        min_entry_period = len(word)
        max_entry_period = 0

        for entry in p:
            _, period, count, _ = entry
            total += period * count

            min_entry_count = min(min_entry_count, count)
            max_entry_count = max(max_entry_count, count)

            min_entry_period = min(min_entry_period, period)
            max_entry_period = max(max_entry_period, period)

        gaps = path2gaps(p, word)
        print(f"length: {len(word)}")
        print(f"cover: {total}")
        print(f"repeat units: {len(p)}")
        print(f"gaps: {gaps}")
        print(f"count: {min_entry_count} - {max_entry_count}")
        print(f"period: {min_entry_period} - {max_entry_period}")

        min_count = min(min_count, min_entry_count)
        max_count = max(max_count, max_entry_count)

        min_period = min(min_period, min_entry_period)
        max_period = max(max_period, max_entry_period)

        print()

    print("Min repeat units:", min_repeats)
    print("Max repeat units:", max_repeats)

    print("Min count:", min_count)
    print("Max count:", max_count)

    print("Min period:", min_period)
    print("Max period:", max_period)


def unique_pmrs(cover):
    pmrs = [x[3] for x in cover]
    if sorted(pmrs) == sorted(set(pmrs)):
        return True
    return False


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
    max_cover = cover(word, pmrs)
    print_tables(n, word, inv, max_cover)

    paths = inv2paths(inv, pmrs, max_cover, overlap, word)
    print([path2hgvs(path, word) for path in paths if cover_length(path) == max_cover[-1]])
    paths = inv2paths(inv, pmrs, max_cover, overlap, word)
    print([path2hgvs(path, word) for path in paths if unique_pmrs(path) and cover_length(path) == max_cover[-1]])
    paths = inv2paths(inv, pmrs, max_cover, overlap, word)
    print(set([path2hgvs(path, word) for path in paths if unique_pmrs(path) and cover_length(path) == max_cover[-1]]))
    # metrics(paths, word, pmrs)


if __name__ == "__main__":
    main()
