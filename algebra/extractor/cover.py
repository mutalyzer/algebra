from os.path import commonprefix
from itertools import combinations, product
import sys


def find_pmrs(word):
    n = len(word)
    pmrs = []
    inv = [None] * n
    for period in range(1, n // 2 + 1):
        start = 0
        while start < n - period:
            pattern = word[start:start + period]

            if inv[start] is not None:
                pmr_start, pmr_period, pmr_count, _ = pmrs[inv[start]]
                if word[pmr_start:pmr_start + pmr_period * pmr_count].startswith(2 * pattern):
                    start += pmr_period * pmr_count - period + 1
                    continue

            count = 1
            while pattern == word[start + period * count:start + period * (count + 1)]:
                count += 1

            remainder = len(commonprefix([pattern, word[start + period * count:]]))
            if count > 1:
                inv[start] = len(pmrs)
                pmrs.append((start, period, count, remainder))
                start += period * (count - 1) + 1
            else:
                start += remainder + 1

    return pmrs


def pmr_interval(pmr):
    start, period, count, remainder = pmr
    return start + 2 * period - 1, start + count * period + remainder


def inv_array(n, pmrs):
    inv = [[] for _ in range(n)]
    for idx, pmr in enumerate(pmrs):
        for pos in range(*pmr_interval(pmr)):
            inv[pos].append(idx)
    return inv


def overlapping(pmrs):
    overlap = [0] * len(pmrs)
    for rhs_idx, rhs in enumerate(pmrs[1:], start=1):
        for lhs in pmrs[:rhs_idx]:
            lhs_start, lhs_period, lhs_count, lhs_remainder = lhs
            rhs_start, rhs_period, rhs_count, rhs_remainder = rhs
            lhs_end = lhs_start + lhs_period * lhs_count + lhs_remainder
            rhs_end = rhs_start + rhs_period * rhs_count + rhs_remainder
            if lhs_start < rhs_start < lhs_end < rhs_end:
                overlap[rhs_idx] = max(overlap[rhs_idx], lhs_end)
    return overlap


def cover(word, pmrs, inv=None, overlap=None):
    n = len(word)
    if not inv:
        inv = inv_array(n, pmrs)
    if not overlap:
        overlap = overlapping(pmrs)

    max_cover = [0] * n
    for pos in range(1, n):
        value = max_cover[pos - 1]

        for idx in inv[pos]:
            start, period, *_ = pmrs[idx]
            count = (pos - start + 1) // period
            length = period * count

            prev_value = 0
            if pos - length > 0:
                prev_value = max_cover[pos - length]

            value = max(value, prev_value + length)

            for p in range(start, min(overlap[idx], pos - period * 2 + 1)):
                count = (pos - p) // period
                length = period * count

                prev_value = max_cover[pos - length]
                value = max(value, prev_value + length)

        max_cover[pos] = value

    return max_cover


def brute_cover(word, pmrs, prev=None, n=0, max_cover=0):
    def intersect(lhs, rhs):
        lhs_start, lhs_period, lhs_count = lhs
        rhs_start, *_ = rhs
        return rhs_start < lhs_start + lhs_period * lhs_count

    if n >= len(pmrs):
        return max_cover

    local = brute_cover(word, pmrs, prev, n + 1, max_cover)

    start, period, count, remainder = pmrs[n]
    for i in range(remainder + 1):
        if prev is None or not intersect(prev, (start + i, period, count)):
            local = max(local, brute_cover(word, pmrs, (start + i, period, count), n + 1, max_cover + period * count))

    for i in range(remainder + period + 1):
        for j in range(2, count):
            for k in range(count - j):
                if prev is None or not intersect(prev, (start + i + period * k, period, j)):
                    local = max(local, brute_cover(word, pmrs, (start + i + period * k, period, j), n + 1, max_cover + period * j))

    return local


def cover_length(cover):
    return sum([period * count for _, period, count, *_ in cover])


def brute_cover_alt(word, pmrs):
    def intersect(a, b):
        a_start, a_period, a_count = a
        b_start, *_ = b
        return b_start < a_start + a_period * a_count

    def bcover(pmrs, n, cover):
        print(f"n: {n}, cover: {cover}")
        if n >= len(pmrs):
            print(f"result: {cover}")
            yield list(cover)
            return

        # without this pmr
        yield from bcover(pmrs, n + 1, cover)

        prev = None
        if len(cover):
            prev = cover[-1]

        start, period, count, remainder = pmrs[n]
        for i in range(remainder + 1):
            # the complete pmr
            entry = start + i, period, count
            print(f"entry 1 {entry} {word[start + i: start + i + period * count]}")
            if prev is None or not intersect(prev, entry):
                cover.append(entry)
                yield from bcover(pmrs, n + 1, cover)
                # yield from bcover(pmrs, n - 1, cover)
                cover.pop()
            else:
                print("skipped because of intersect")

        print("yo")

        for i in range(remainder + period + 1):
            for j in range(2, count):
                for k in range(count - j):
                    entry = start + i + period * k, period, j
                    print(f"entry 2 {entry}")
                    if prev is None or not intersect(prev, entry):
                        cover.append(entry)
                        yield from bcover(pmrs, n + 1, cover)
                        cover.pop()
                    else:
                        print("skipped because of intersect")

    return bcover(pmrs, 0, [])


def cartesian_cover(pmrs):
    def intersect(lhs, rhs):
        lhs_start, lhs_period, lhs_count, lhs_pmrs = lhs
        lhs_end = lhs_start + lhs_period * lhs_count - 1
        rhs_start, rhs_period, rhs_count, rhs_pmrs = rhs

        if lhs_start < rhs_start:
            return lhs_end >= rhs_start

        rhs_end = rhs_start + rhs_period * rhs_count - 1
        return rhs_end >= lhs_start

    def complete(pmrs, n):
        # print(f"pmrs: {n}")
        start, period, count, remainder = pmrs[n]

        yield
        for offset in range(period):
            # print(f"offset: {offset}")

            count2 = count
            if offset > remainder:
                count2 -= 1

            for prefix in range(count2 - 1):
                # print(f"prefix: {prefix}")

                for c in range(2, count2 - prefix + 1):
                    # print(f"count: {c}")

                    entry = start + offset + prefix * period, period, c, n
                    # print(f"entry:", entry)
                    yield entry

    lists = [complete(pmrs, n) for n in range(len(pmrs))]
    for candidate in product(*lists):
        filtered = list(filter(None, candidate))
        # print("candidate:", filtered)

        if len(filtered) == 0:
            continue

        if len(filtered) == 1:
            # print("single solution:", list(filtered))
            yield sorted(filtered)
            continue

        for a, b in combinations(filtered, 2):
            # print("combo:", a, b)
            if intersect(a, b):
                # print("intersect, break")
                break
            # print("No intersect")
        else:
            # print("solution:", filtered)
            yield sorted(filtered)


def print_array(array):
    for element in array:
        print(f"{element:3}", end="")
    print()


def print_tables(n, word, inv, max_cover):
    print_array(range(n))
    print("  ", end="")
    for ch in word:
        print(f"{ch:3}", end="")
    print()
    for y in range(max([len(x) for x in inv])):
        for x in inv:
            if len(x) <= y:
                if y == 0:
                    print("  .", end="")
                else:
                    print("   ", end="")
            else:
                print(f"{x[y]:3}", end="")
        print()
    print_array(max_cover)


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
    print(brute_cover(word, pmrs))

    # covers = list(brute_cover_alt(word, pmrs))
    # bmax = max(map(cover_length, covers))
    from walker import path2hgvs
    # print({path2hgvs(c, word) for c in covers if cover_length(c) == bmax})

    # print(list(cartesian_cover(word, pmrs, max_cover[-1])))
    covers = list(cartesian_cover(pmrs))
    bmax = max(map(cover_length, covers))
    # print({path2hgvs(c, word) for c in covers})
    print({path2hgvs(c, word) for c in covers if cover_length(c) == bmax})


if __name__ == "__main__":
    main()
