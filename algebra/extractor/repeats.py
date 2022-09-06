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


def subswithseqs(subs, word):
    for start, period, count, remainder in subs:
        print(str((start, period, count, remainder)) + ",  #", word[start:start + period], word[start:start + period] * count)


def main():
    if len(sys.argv) < 2:
        return -1
    input = sys.argv[1]

    # Get all the repeats as (position, period, count, remainder) quadruples
    subs = extract_repeats(input)
    print(subs)
    print(len(subs))
    subswithseqs(subs, input)


if __name__ == '__main__':
    main()
