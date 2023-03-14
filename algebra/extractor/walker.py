import sys
from algebra.extractor.cover import cover, cover_length, find_pmrs, inv_array, overlapping, print_tables


def kanji(max_cover, path, inv, pmrs, word):
    elements, total = kfill(max_cover, sorted(path), inv, pmrs, word)
    return "".join(elements), total


def kfill(max_cover, subs, inv, pmrs, word):
    desc = []

    end = -1
    total = 0
    for start, period, count, *_ in subs:

        # This position already covered by previous iteration
        if end >= start:
            continue

        print(f"start: {start}, period: {period}, count: {count}")

        # Add leading bases that were not part of a repeat
        if end + 1 < start:
            desc.append(word[end + 1: start])
            total += start - end - 1
            print(f"add {start - end - 1} in the begin")

        # Add repeat to description
        if period * count > 2:
            desc.append(f"({word[start:start + period]})^{count}")
            total += period + 1
            print(f"add {period + 1} because of repeat")
        else:
            desc.append(word[start:start + period * count])
            total += period * count
            print(f"add {period * count} because of small repeat")

        # Set new end base
        end = start + period * count - 1

    # Add trailing bases that was not part of a repeat
    if end + 1 < len(word):
        desc.append(word[end + 1:])
        total += len(word) - end - 1
        print(f"add {len(word) - end - 1} at end")

    return desc, total


def walk(inv, pmrs, max_cover, pos, path, start):
    # Done!
    if pos < start:
        yield path
        return

    # Move to the next position if there is an (equally) maximal solution
    if max_cover[pos - 1] == max_cover[pos]:
        yield from walk(inv, pmrs, max_cover, pos - 1, path, start)

    # Try all pmrs at this position
    for idx in inv[pos]:
        if idx in [x[3] for x in path]:
            # Skip if we already used this pmrs for this path
            continue
            # pass

        begin, period, _, _ = pmrs[idx]

        # Don't let candidate collide with previous entry
        max_count = (pos - begin + 1) // period

        # Try all counts downwards
        for count in range(max_count, 1, -1):
            # Move left until we find a position with a solution
            prev_pos = pos - period * count
            while not inv[prev_pos]:
                prev_pos -= 1
                if prev_pos < 1:
                    break

            # Move to previous position if there is an (equally) maximal solution (or if we are the last)
            if max_cover[max(0, prev_pos)] + period * count == max_cover[pos]:
                # The actual entry
                entry = pos - period * count + 1, period, count, idx
                yield from walk(inv, pmrs, max_cover, prev_pos, path + [entry], start)


def walk2(inv, pmrs, max_cover, pos, path, start, word):
    # Done!
    if pos <= start:
        print(f"done: {path}")
        if pos == start:
            yield word[pos] + path
        else:
            yield path
        return


    # Move to the next position if there is an (equally) maximal solution
    if max_cover[pos - 1] == max_cover[pos]:
        print(f"word pos: {word[pos]}")
        yield from walk2(inv, pmrs, max_cover, pos - 1, word[pos] + path, start, word)

    # Try all pmrs at this position
    for idx in inv[pos]:
        print(f"loop {idx}")
        # if idx in [x[3] for x in path]:
        #     # Skip if we already used this pmrs for this path
        #     continue
        #     # pass

        begin, period, _, _ = pmrs[idx]

        # Don't let candidate collide with previous entry
        max_count = (pos - begin + 1) // period

        # Try all counts downwards
        for count in range(max_count, 1, -1):
            # Move left until we find a position with a solution
            prev_pos = pos - period * count
            while not inv[prev_pos]:
                prev_pos -= 1
                # if prev_pos < 0:
                #     break

            # Move to previous position if there is an (equally) maximal solution (or if we are the last)
            if max_cover[max(0, prev_pos)] + period * count == max_cover[pos]:
                # The actual entry
                entry = pos - period * count + 1, period, count, idx
                print(f"entry: {entry}")
                for entry in walk2(inv, pmrs, max_cover, begin + period - 1, '', begin, word):
                    print("recurse")
                    yield from walk2(inv, pmrs, max_cover, prev_pos, f"({entry})^{count}" + path, start, word)


class Node:
    def __init__(self, start, period, count, pmrs):
        self.start = start
        self.period = period
        self.count = count
        self.pmrs = pmrs

        self.children = []

    def __str__(self):
        return f"{self.start}, {self.period}, {self.count}, {self.pmrs}"


def walk3(inv, pmrs, max_cover, pos, node, start, word):
    # Done!
    if pos <= start:
        print(f"done: {node} ({pos})")
        yield node
        return

    # Move to the next position if there is an (equally) maximal solution
    if max_cover[pos - 1] == max_cover[pos]:
        yield from walk3(inv, pmrs, max_cover, pos - 1, node, start, word)

    # Try all pmrs at this position
    for idx in inv[pos]:
        print(f"loop {idx} pos: {pos}")
        # if idx in [x[3] for x in path]:
        #     # Skip if we already used this pmrs for this path
        #     continue
        #     # pass

        begin, period, _, _ = pmrs[idx]

        # Don't let candidate collide with previous entry
        max_count = (pos - begin + 1) // period

        # Try all counts downwards
        for count in range(max_count, 1, -1):
            # Move left until we find a position with a solution
            prev_pos = pos - period * count
            while not inv[prev_pos]:
                prev_pos -= 1
                if prev_pos < 0:
                    break

            # Move to previous position if there is an (equally) maximal solution (or if we are the last)
            if max_cover[max(0, prev_pos)] + period * count == max_cover[pos]:
                # The actual entry
                entry = Node(pos - period * count + 1, period, count, idx)
                print(f"entry: {entry}")
                print("going to recurse", entry)
                for x in walk3(inv, pmrs, max_cover, begin + period - 1, entry, begin, word):
                    print(f"x: {x}")
                    child = walk3(inv, pmrs, max_cover, prev_pos, x, start, word)
                    print(f"child: {list(child)[0]}")
                    print("Adding child")
                    # entry.children.append(child)


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
    return ";".join(fill(sorted(path), word)[0])


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
    indices = [e[3] for e in cover]
    if len(indices) == len(set(indices)):
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
    print(inv)
    max_cover = cover(word, pmrs)
    print(max_cover)
    print_tables(n, word, inv, max_cover)

    paths = walk(inv, pmrs, max_cover, len(max_cover) - 1, [], 1)
    for path in paths:
        desc, total = kanji(max_cover, path, inv, pmrs, word)
        print(total, desc)

    # paths = walk3(inv, pmrs, max_cover, len(max_cover) - 1, None, 1, word)
    #print([path2hgvs(path, word) for path in paths])
    # print(list(paths))


    # paths = list(walk(inv, pmrs, max_cover, overlap, len(max_cover) - 1, []))
    # print([path for path in paths if cover_length(path) == max_cover[-1]])
    # print([path2hgvs(path, word) for path in paths if cover_length(path) == max_cover[-1]])
    # paths = walk(inv, pmrs, overlap, len(max_cover) - 1, [])
    # print([path2hgvs(path, word) for path in paths if unique_pmrs(path) and cover_length(path) == max_cover[-1]])
    # paths = walk(inv, pmrs, overlap, len(max_cover) - 1, [])
    # print(set([path2hgvs(path, word) for path in paths if unique_pmrs(path) and cover_length(path) == max_cover[-1]]))
    # paths = walk(inv, pmrs, overlap, len(max_cover) - 1, [])
    # metrics(paths, word, pmrs)
    # paths = walk(inv, pmrs, overlap, len(max_cover) - 1, [])
    # metrics([path for path in paths if unique_pmrs(path) and cover_length(path) == max_cover[-1]], word, pmrs)


if __name__ == "__main__":
    main()
