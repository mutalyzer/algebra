from os.path import commonprefix
import sys


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

    pred = -1
    for pos in range(1, n):
        # default to the previous value ("O-class")
        value = max_cover[pos - 1]
        # print("pos", pos)

        # if len(inv[pos]) == 0:
        #     pred = -1

        for idx in inv[pos]:
            start, period, count, remainder = pmrs[idx]

            overlap = 0
            if pred != -1:
                pred_start, pred_period, pred_count, pred_remainder = pmrs[pred]
                pred_end = pred_start + pred_period * pred_count + pred_remainder
                end = start + period * count + remainder

                if pred_end > start and pred_end < end and pred_start < start:
                    overlap = pred_end - start
                    print(f"{idx} overlaps {pred} at {pos}: {overlap}")
                else:
                    pred = idx
            else:
                pred = idx

            # "N-class"
            # print("N class")
            real_count = (pos - start + 1) // period
            real_length = real_count * period

            prev_value = 0
            if pos - real_length >= 0:  # > 0 ?
                prev_value = max_cover[pos - real_length]

            value = max(value, prev_value + real_length)
            # print("value", value)

            # "N-1-class"
            # print("N-1 class")
            real_count -= overlap // period

            if real_count > 1:
                real_length = real_count * period
                # print("real", real_count, real_length)

                prev_value = 0
                if pos - real_length >= 0:
                    prev_value = max_cover[pos - real_length]

                # print("prev", pos - real_length, prev_value)

                value = max(value, prev_value + real_length)
                # print("value", value)


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


def main():
    if len(sys.argv) < 2:
        return
    word = sys.argv[1]
    n = len(word)
    print(n, word)
    pmrs = find_pmrs(word)
    print("pmrs", pmrs)
    inv = inv_array(n, pmrs)
    print("inv", inv)
    max_cover = cover(word, pmrs)
    print("max_cover", max_cover)
    #for e in extract_cover(len(word) - 1, pmrs, inv, max_cover, []):
    #    print(to_hgvs(word, e))

    #print("start bcover")
    #print(sorted([(start, period, count, remainder) for start, period, count, remainder in pmrs]))
    #print(brute_cover(word, sorted([(start, period, count, remainder) for start, period, count, remainder in pmrs])))
    #print("end bcover")


if __name__ == "__main__":
    main()
