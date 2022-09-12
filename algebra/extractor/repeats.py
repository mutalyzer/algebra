from itertools import chain, combinations
from os.path import commonprefix
import sys


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


def subswithseqs(subs, word):
    for start, period, count, remainder in subs:
        print(str((start, period, count, remainder)) + ",  #", word[start:start + period], word[start:start + period] * count)


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


def brutepower(pmrs, word):

    def intersect(a, b):
        a_start, a_period, a_count, a_remainder = a
        a_end = a_start + a_period * a_count

        b_start, b_period, b_count, b_remainder = b

        return b_start < a_end

    pmrs_complete = []
    for start, period, count, remainder in pmrs:
        pmrs_complete.append((start, period, count, 0))  # add self
        length = period * count + remainder
        for l in range(length - 1, 1, -1):
            for offset in range(length - l + 1):
                if l // period > 1:
                    pmrs_complete.append((start + offset, period, l // period, 0))

    # print(pmrs_complete)

    y = 0
    solutions = []
    for subset in powerset(pmrs_complete):
        # print(subset)
        l = list(subset)


        conflict = False
        for idx in range(1, len(l)):
            if intersect(l[idx - 1], l[idx]):
                # print("conflict", l[idx - 1], l[idx])
                conflict = True
                break

        if not conflict:
            z = 0
            for x in l:
                if x[2] > 1:
                    z += x[1] * x[2]
            # print("answer:", l, z)
            old_max = y
            y = max(y, z)

            if z == y:
                if z > old_max:
                    solutions = [l]
                else:
                    solutions.append(l)

    print("solutions", solutions)

    return y, solutions


def main():
    if len(sys.argv) < 2:
        return -1
    word = sys.argv[1]

    # Get all the repeats as (position, period, count, remainder) quadruples
    runs = extract_repeats(word)
    print(runs)
    print(len(runs))
    subswithseqs(runs, word)
    _, solutions = brutepower(runs, word)
    for sol in solutions:
        print("Solution")
        subswithseqs(sol, word)


if __name__ == '__main__':
    main()
