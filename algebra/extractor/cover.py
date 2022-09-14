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

    values = [0] * len(pmrs)
    ends = [0] * len(pmrs)
    max_cover = [0] * n

    hwm = 0
    for pos in range(1, n):
        # default to the previous value
        value = max_cover[pos - 1]
        print(f"pos: {pos}")

        for idx in inv[pos]:
            print(f"idx: {idx}")
            start, period, count, remainder = pmrs[idx]

            if hwm >= start:
                real_count = (pos - hwm) // period
                if real_count > 1:
                    real_length = real_count * period
                    if max_cover[hwm] + real_length > value:
                        value = max_cover[hwm] + real_length
                        if value > values[idx] and (values[idx] == 0 or len(inv[pos]) == 1):
                            values[idx] = value
                            ends[idx] = pos
                            print(f"Update because of hwm values[idx] {values[idx]}")
                            print(f"Update because of hwm ends[idx] {ends[idx]}")

            real_count = (pos - start + 1) // period
            if real_count > 1:
                real_length = real_count * period
                prev_value = 0
                if pos - real_length >= 0:
                    prev_value = max_cover[pos - real_length]

                if prev_value + real_length > value:
                    value = prev_value + real_length
                    if value > values[idx] and (values[idx] == 0 or len(inv[pos]) == 1):
                        values[idx] = value
                        ends[idx] = pos
                        print(f"Update values[idx] {values[idx]}")
                        print(f"Update ends[idx] {ends[idx]}")

            if start + period * count + remainder - 1 == pos and values[idx] > 0:
                print("update hwm (pos, idx, ends[idx])", pos, idx, ends[idx])
                hwm = ends[idx]

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


def main():
    word = fib_word(9)
    n = len(word)
    pmrs = find_pmrs(word)
    print(len(pmrs))
    for idx, pmr in enumerate(pmrs):
        print(f"        {pmr},  # {idx:2}: {word[pmr[0]:pmr[0] + pmr[1]]}")


def to_hgvs(word, repeats):
    def hgvs():
        start = 0
        for repeat in repeats:
            if repeat[0] > start:
                yield word[start:repeat[0]]
            yield f"{word[repeat[0]:repeat[0] + repeat[1]]}[{repeat[2]}]"
            start = repeat[0] + repeat[1] * repeat[2]
        if start < len(word):
            yield word[start:]
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


def main():
    # word = fib_word(9)
    # n = len(word)
    # pmrs = find_pmrs(word)
    # print(len(pmrs))
    # for idx, pmr in enumerate(pmrs):
    #     print(f"        {pmr},  # {idx:2}: {word[pmr[0]:pmr[0] + pmr[1]]}")
    #
    # inv = inv_array(n, pmrs)
    # print(n, inv, sum([len(i) for i in inv]))
    # print(cover(word, pmrs))
    # return

    #print(find_pmrs("CCCCACCAT"))
    #return

    for word, pmrs, inv, max_cover, hgvs in TESTS[4:5]:
        print(word)
        # my_pmrs = find_pmrs(word)
        # assert set(pmrs) == set(my_pmrs)
        # continue

        n = len(word)
        # print("AAA", cover(word, pmrs))
        assert inv == inv_array(n, pmrs)
        assert max_cover == cover(word, pmrs)

        print()
        print()
        print()
        for e in extract_cover(len(word) - 1, pmrs, inv, max_cover, []):
            print(to_hgvs(word, e))


if __name__ == "__main__":
    main()
