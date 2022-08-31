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
                    # Don't continue for less/contained repeats.
                    break

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


def extract_global(subs, word):
    desc = []

    end = -1
    window_end = -1
    idx2 = -1

    for idx, (start, period, count) in enumerate(subs):
        print(f"idx: {idx} idx2: {idx2} window_end: {window_end}")
        print(f"idx: {idx} start: {start} period: {period} count: {count} span: {period * count}")
        max_span = period * count
        max_idx = idx

        if idx < idx2:
            print("Continue because of index")
            continue

        if start <= window_end:
            print("Continue because of window_end")
            continue

        window_end = start + period * count - 1
        print(f"Window end: {window_end} span: {period * count}")

        idx2 = idx + 1

        stack = [idx]

        next_end = window_end
        while idx2 < len(subs):
            next_start, next_period, next_count = subs[idx2]
            next_span = next_period * next_count
            print(f"idx2: {idx2} start: {next_start} period: {next_period} count: {next_count} span: {next_span}")

            if next_start <= window_end:
                print(f"Overlap. Span: {next_span}")
                if next_span > max_span:
                    print("Higher span!")
                    max_idx = idx2
                    max_span = next_span
                    next_end = next_start + next_span - 1
                    stack.append(idx2)
                    print(f"stack: {stack}")
                elif next_span < max_span:
                    print("Lower span!")
                else:
                    print("Same span!")
            else:
                print("No overlap")
                window_end = next_end
                break

            print(f"max_idx: {max_idx} max_span: {max_span} (end of while)")
            # max_start, max_period, max_count = subs[max_idx]
            # desc.append(f"{word[max_start:max_start + max_period]}[{max_count}]")

            idx2 += 1
        print(f"stack: {stack} after while")

        # for stack_idx in reversed(stack):
        #     stack_start, stack_period, stack_count = subs[stack_idx]
        #     print(stack_start, stack_period, stack_count)
        # max_start, max_period, max_count = subs[max_idx]
        # desc.append(f"{word[max_start:max_start + max_period]}[{max_count}]")
        stack_start, stack_period, stack_count = subs[stack[-1]]
        print(stack_start, stack_period, stack_count)
        desc.append(f"{word[stack_start:stack_start + stack_period]}[{stack_count}]")

        print()

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

    print()

    # Extract HGVS-like descriptions
    desc = extract_global(subs, input)
    print(";".join(desc))

if __name__ == '__main__':
    main()
