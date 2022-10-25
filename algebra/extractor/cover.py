from os.path import commonprefix
import sys
from itertools import permutations


def pmr_interval(pmr):
    start, period, count, remainder = pmr
    return start + 2 * period - 1, start + count * period + remainder


def inv_array(n, pmrs):
    inv = [[] for _ in range(n)]
    for idx, pmr in enumerate(pmrs):
        for pos in range(*pmr_interval(pmr)):
            inv[pos].append(idx)
    return inv



def cover(word, pmrs, inv=None):
    n = len(word)
    if not inv:
        inv = inv_array(n, pmrs)

    max_cover = [0] * n

    overlap_array = [0] * len(pmrs)
    conflict_array = [-1] * len(pmrs)
    for o, p in permutations(range(len(pmrs)), 2):
        start, period, count, remainder = pmrs[o]
        pred_start, pred_period, pred_count, pred_remainder = pmrs[p]
        end = start + period * count + remainder
        pred_end = pred_start + pred_period * pred_count + pred_remainder

        if pred_start < start < pred_end < end:
            overlap = pred_end - start
            print(f"{o} overlaps {p}: {overlap}")
            if overlap > overlap_array[o]:
                overlap_array[o] = overlap
                conflict_array[o] = p

    print("Conflict array", conflict_array)

    # pred = -1
    for pos in range(1, n):
        # default to the previous value ("O-class")
        value = max_cover[pos - 1]
        print("pos", pos)

        # if len(inv[pos]) == 0:
        #     pred = -1

        for idx in inv[pos]:
            start, period, count, remainder = pmrs[idx]
            print("idx ", idx)
            print("overlap", overlap_array[idx])
            print("conflict", conflict_array[idx])

            # overlap = 0
            # if pred != -1:
            #     pred_start, pred_period, pred_count, pred_remainder = pmrs[pred]
            #     pred_end = pred_start + pred_period * pred_count + pred_remainder
            #     end = start + period * count + remainder
            #
            #     if pred_end > start and pred_end < end and pred_start < start:
            #         overlap = pred_end - start
            #         print(f"{idx} overlaps {pred} at {pos}: {overlap}")
            #     else:
            #         pred = idx
            # else:
            #     pred = idx

            # "N-class"
            print("N class")
            real_count = (pos - start + 1) // period
            real_length = real_count * period

            prev_value = 0
            if pos - real_length >= 0:  # > 0 ?
                prev_value = max_cover[pos - real_length]

            print("value", prev_value + real_length)

            value = max(value, prev_value + real_length)

            if conflict_array[idx] != -1:
                print("Conflict!")
                x_start, x_period, x_count, x_remainder = pmrs[conflict_array[idx]]

                for x in range(start, start + overlap_array[idx]):
                    print("x", x)
                    if x < start:
                        print("skip large x")
                        continue
                    real_count = (pos - x) // period
                    print("real count", real_count)

                    if real_count > 1:
                        real_length = real_count * period
                        print("real length", real_length)

                        prev_value = 0
                        if pos - real_length >= 0:
                            prev_value = max_cover[pos - real_length]

                        value = max(value, prev_value + real_length)
                        print("value", value)

            # "N-1-class"
            # print("N-1 class")
            # real_count -= overlap // period
            #
            # if real_count > 1:
            #     real_length = real_count * period
            #     # print("real", real_count, real_length)
            #
            #     prev_value = 0
            #     if pos - real_length >= 0:
            #         prev_value = max_cover[pos - real_length]
            #
            #     # print("prev", pos - real_length, prev_value)
            #
            #     value = max(value, prev_value + real_length)
            #     # print("value", value)


        max_cover[pos] = value

    return max_cover


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


def fib_word(n):
    fw0 = "A"
    if n <= 0:
        return fw0

    fw1 = "C"
    for _ in range(1, n):
        fw0, fw1 = fw1, fw1 + fw0
    return fw1


def frac(n):
    fw0 = "A"
    if n <= 0:
        return fw0
    fw1 = "AC"
    if n == 1:
        return fw1
    fw2 = "AG"
    if n == 2:
        return fw2
    fw3 = "AT"
    if n == 3:
        return fw3

    for _ in range(3, n):
        fw0, fw1, fw2, fw3 = fw1, fw2, fw3, fw3 + fw2 + fw1 + fw0
    return fw3


def to_hgvs(word, repeats):
    def hgvs():
        pos = 0
        for start, period, count in repeats:
            if start > pos:
                yield word[pos:start]
            yield f"{word[start:start + period]}[{count}]"
            pos = start + period * count
        if pos < len(word):
            yield word[pos:]
    return ";".join(hgvs())


def extract_cover(pos, pmrs, inv, max_cover, hwm, start=0, cover=[]):
    if pos < start:
        print("end pos", pos, cover)
        yield list(cover)
        return

    if not inv[pos]:
        print("inv pos", pos)
        yield from extract_cover(pos - 1, pmrs, inv, max_cover, hwm, start, cover)
        return

    print("work")

    if max_cover[pos] == max_cover[pos - 1]:
        print("yield left", pos)
        yield from extract_cover(pos - 1, pmrs, inv, max_cover, hwm, start, cover)

    for idx in inv[pos]:
        start, period, count, remainder = pmrs[idx]
        length = period * count
        if pos - length == -1 or max_cover[pos] == max_cover[pos - length] + length:
            print("pmrs", pmrs[idx], pos)
            offset = pos + 1 - period * count
            cover.append((start + offset, period, count))
            yield from extract_cover(pos - length, pmrs, inv, max_cover, hwm, start, cover)
            cover.pop()


def brute_cover(word, pmrs):
    def intersect(a, b):
        a_start, a_period, a_count = a
        b_start, *_ = b
        return b_start < a_start + a_period * a_count

    def cover_length(cover):
        return sum([period * count for _, period, count in cover])

    def bcover(pmrs, n=0, cover=[]):
        if n >= len(pmrs):
            yield cover_length(cover)
            return

        # whitout this prm
        yield from bcover(pmrs, n + 1, cover)

        prev = None
        if len(cover):
            prev = cover[-1]

        start, period, count, remainder = pmrs[n]
        for i in range(remainder + 1):
            # the complete prm
            if prev is None or not intersect(prev, (start + i, period, count)):
                cover.append((start + i, period, count))
                yield from bcover(pmrs, n + 1, cover)
                cover.pop()

        for i in range(remainder + period + 1):
            for j in range(2, count):
                for k in range(count - j):
                    if prev is None or not intersect(prev, (start + i + k * period, period, j)):
                        cover.append((start + i + k * period, period, j))
                        yield from bcover(pmrs, n + 1, cover)
                        cover.pop()

    return max(bcover(pmrs))


def overlapping(pmrs):
    def overlap(pmr_lhs, pmr_rhs):
        start_lhs, period_lhs, count_lhs, remainder_lhs = pmr_lhs
        end_lhs = start_lhs + period_lhs * count_lhs + remainder_lhs
        start_rhs, period_rhs, count_rhs, remainder_rhs = pmr_rhs
        end_rhs = start_rhs + period_rhs * count_rhs + remainder_rhs

        if start_lhs < start_rhs < end_lhs < end_rhs:
            return end_lhs - start_rhs
        return 0

    result = [[] for _ in range(len(pmrs))]
    for idx_rhs, rhs in enumerate(pmrs[1:], start=1):
        if rhs[0] == pmrs[0][0]:
            continue
        for idx_lhs, lhs in enumerate(pmrs[:idx_rhs]):
            length = overlap(lhs, rhs)
            if length:
                result[idx_rhs].append((idx_lhs, length))
    return result


def print_tables(n, word, inv, cover):
    for pos in range(n):
        print(f"{pos:3}", end="")
    print("\n", end="  ")
    for ch in word:
        print(f"{ch:3}", end="")
    print()
    for y in range(max([len(x) for x in inv])):
        for x in inv:
            if len(x) <= y:
                if y == 0:
                    print(f"  .", end="")
                else:
                    print("   ", end="")
            else:
                print(f"{x[y]:3}", end="")
        print()
    for val in cover:
        print(f"{val:3}", end="")
    print()


def main():
    if len(sys.argv) < 2:
        print(f"usage: {sys.argv[0]} word", file=sys.stderr)
        return

    word = sys.argv[1]
    n = len(word)
    print(n, word)

    pmrs = sorted(find_pmrs(word))
    for idx, pmr in enumerate(pmrs):
        print(f"{pmr},  # {idx:2}: {word[pmr[0]:pmr[0] + pmr[1]]}")

    print(overlapping(pmrs))

    inv = inv_array(n, pmrs)
    max_cover = cover(word, pmrs)

    print_tables(n, word, inv, max_cover)
    print(brute_cover(word, pmrs))


if __name__ == "__main__":
    main()
