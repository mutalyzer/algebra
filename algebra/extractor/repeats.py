import sys


def repeats(word):
    results = []

    # For all bases as start position, except the last
    for start in range(len(word)):

        # For all the periods remaining
        for period in range(1, len(word) - start + 1):

            # For all extensions to the right remaining in reverse order
            for extension in reversed(range(period, len(word) - start - period + 1, period)):

                # Determine the number of (potential) repeats
                count = extension // period + 1

                # Add to the results if the haystack is made out of needles exclusively
                if word[start:start + period] * count == word[start:start + period + extension]:
                    results.append((start, period, count))

    return results


def extract(subs, word):
    desc = []

    end = -1
    for start, period, count in subs:

        # This position already covered by previous iteration
        if end >= start:
            continue

        # Add leading bases that were not part of a repeat
        if end + 1 < start:
            desc.append(word[end + 1: start])

        # Add repeat to description
        desc.append(word[start:start + period] + f"[{count}]")

        # Set new end base
        end = start + period * count - 1

    # Add trailing bases that was not part of a repeat
    if end + 1 < len(word):
        desc.append(word[end + 1:])

    return desc


def main():
    if len(sys.argv) < 2:
        return -1
    input = sys.argv[1]

    # Get all the repeats as (position, period, count) triples
    subs = repeats(input)
    print(subs)

    # Extract HGVS-like descriptions
    desc = extract(subs, input)
    print(";".join(desc))


if __name__ == '__main__':
    main()
