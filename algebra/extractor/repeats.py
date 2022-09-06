import sys
from os.path import commonprefix


def extract_repeats(word):
    results = []

    # For all bases as start position, except the last
    for start in range(len(word)):

        # Deal with rotations that belong to us
        if results:
            last_start, last_period, last_count, last_remainder = results[-1]
            if last_remainder > 0 and start == last_start + 1:
                continue

        # For all the periods remaining
        for period in range(1, len(word) - start + 1):

            # For all extensions to the right remaining in reverse order
            for extension in reversed(range(period, len(word) - start - period + 1, period)):

                # Determine the number of (potential) repeats
                count = extension // period + 1

                # Add to the results if the haystack is made out of needles exclusively
                if word[start:start + period] * count == word[start:start + period + extension]:
                    remainder = len(commonprefix([word[start:start + period], word[start + period * count: start + period * (count + 1)]]))

                    back = 1
                    skip = False
                    while len(results) >= back:
                        last_start, last_period, last_count, last_remainder = results[-back]
                        back += 1
                        last_end = last_start + last_period * last_count
                        end = start + period * count
                        if last_end + last_remainder == end + remainder and period % last_period == 0:
                            # Non-primitive, same-period or sub-pattern
                            skip = True
                            break

                    if not skip:
                        results.append((start, period, count, remainder))

                    # Don't continue for less/contained repeats.
                    break

    return results


def koffer(pmrs, word):

    inv = [[] for _ in range(len(word))]

    for idx, (start, period, count, remainder) in enumerate(pmrs):
        begin = start + period * 2 - 1
        end = start + period * count + remainder

        print(idx, begin, end)

        for pos in range(begin, end):
            inv[pos].append(idx)

    for pos, a in enumerate(inv):
        print(pos, a)

    max_cover = [0] * len(word)

    for idx in range(len(max_cover)):
        print("idx", idx)
        prev = 0
        if idx > 0:
            prev = max_cover[idx - 1]

        if not inv[idx]:
            max_cover[idx] = prev
            continue

        # I / 0 = prev
        m = prev
        print("default m + 0", m)
        for pmr in inv[idx]:
            start, period, count, remainder = pmrs[pmr]
            length = period * ((idx - start + 1) // period)
            print("length", length)

            # II / -1
            length2 = length - period
            print("length II", length2)
            count2 = length2 // period
            print("count II", count2)
            idx2 = idx - length2
            print("index II", idx2)
            if count2 > 1:
                m = max(m, max_cover[idx2] + length2)
            else:
                m = max(m, max_cover[idx2])
            print("m after II", m)

            # III / F
            if idx - length > 0:
                m = max(m, max_cover[idx - length] + length)
            else:
                m = max(m, length)

        print("m", m)
        max_cover[idx] = m

    for pos, a in enumerate(max_cover):
        print(pos, a)

    return max_cover[-1]


def subswithseqs(subs, word):
    for start, period, count, remainder in subs:
        print(str((start, period, count, remainder)) + ",  #", word[start:start + period], word[start:start + period] * count)


def main():
    if len(sys.argv) < 2:
        return -1
    word = sys.argv[1]

    # Get all the repeats as (position, period, count, remainder) quadruples
    runs = extract_repeats(word)
    print(runs)
    print(len(runs))
    subswithseqs(runs, word)
    koffer(runs, word)


if __name__ == '__main__':
    main()
