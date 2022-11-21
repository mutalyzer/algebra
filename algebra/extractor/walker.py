import sys
from algebra.extractor.cover import brute_cover, cover, find_pmrs, inv_array, overlapping, print_tables


def walk(solutions, max_result, pos=None, path=None):
    # print(f"walk: {pos}")
    for idx, period, cnt in solutions[pos]:
        for count in range(cnt, 1, -1):
            # print(f"count: {count}")
            # Adjacent entries from same pmrs are not part of minimal solution
            if path and idx == path[-1][3]:
                # print(f"same pmrs: {idx}")
                break  # TODO: continue?!

            # The actual entry
            entry = pos - period * count + 1, period, count, idx
            # print(f"entry: {entry}")

            # Move left until we find a position with a solution
            prev_pos = pos - period * count
            while not solutions[prev_pos]:
                # print("lowering prev_pos")
                prev_pos -= 1
                if prev_pos < 1:
                    break

            # There is possibly more to explore
            if prev_pos > 0:
                yield from walk(solutions, max_result, prev_pos, path + [entry])
            # Otherwise, yield result
            else:
                # print(f"result: {path + [entry]}")
                if sum([x[1] * x[2] for x in path + [entry]]) == max_result:
                    yield reversed(path + [entry])


def walk3(inv, pmrs, max_cover, overlap, word, pos, path):
    print("walk3")

    print("pos", pos)
    for idx in inv[pos]:
        begin = pmrs[idx][0]
        print(f"begin: {begin}")
        period = pmrs[idx][1]
        print(f"period: {period}")
        cnt = pmrs[idx][2]
        print(f"cnt: {cnt}")

        if path and idx == path[-1][3]:
            # print(f"same pmrs: {idx}")
            continue  # TODO: continue?!

        for count in range(cnt, 1, -1):
            print(f"count: {count}")

            start = max(pos - period * count + 1, begin)
            print(f"start: {start}")

            if path and start + count * period > path[-1][0]:
                print("Skip, doesn't fit")
                continue

            entry = start, period, count, idx
            print(f"entry: {entry}")

            # Move left until we find a position with a solution
            prev_pos = pos - period * count
            while not inv[prev_pos]:
                print("lowering prev_pos")
                prev_pos -= 1
                if prev_pos < 1:
                    break

            # There is possibly more to explore
            if prev_pos > 0:
                yield from walk3(inv, pmrs, max_cover, overlap, word, prev_pos, path + [entry])
            # Otherwise, yield result
            else:
                if sum([x[1] * x[2] for x in path + [entry]]) == max_cover[-1]:
                    print(f"result: {path + [entry]}")
                    yield reversed(path + [entry])
                else:
                    print("Skipping non-maximal solution")


def walk2(inv, pmrs, max_cover, overlap, word, pos, result):
    print("walk2")
    if pos < 0:
        if sum([x[1] * x[2] for x in result]) < max_cover[-1]:
            print("Non-maximal end result")
            pass
            # yield reversed(result)
        else:
            print(f"done!! at {pos}: {result}")
            yield reversed(result)

        return
    print(f"result: {result}")

    while not inv[pos]:
        print("Decreasing pos")
        pos -= 1

    # if max_cover[pos] == max_cover[pos - 1]:
    #     print(f"Position {pos} not part of solution")
    #     return

    print("pos", pos)
    for idx in inv[pos]:
        print("pos", pos)
        print("pmrs", idx, pmrs[idx])
        print("overlap", overlap[idx])

        if result:
            prev_pos, prev_period, prev_count, prev_idx = result[-1]
            if prev_idx == idx:
                print("No entry because of same pmrs")
                yield from walk2(inv, pmrs, max_cover, overlap, word, pos - 1, result)
                # continue

        begin = pmrs[idx][0]
        length = pmrs[idx][1]
        if not overlap[idx]:
            end = pos - 2 * length + 1
        else:
            end = overlap[idx]
            print(f"end: {end} because of overlap")

        print(f"begin: {begin} end: {end} length: {length}")
        for start in range(end, begin - 1, -length):
            print(f"X start: {start}")
            value = pos - start + 1  # % (pos - start + 1) % length
            print("value", value)
            count = value // length

            # if start > 0:
            #     prev = max_cover[start - 1]
            # else:
            #     prev = 0
            prev = max_cover[max(start - 1, 0)]
            print(f"prev: {prev}")

            if value + prev != max_cover[pos]:
                print("No max value")
                continue

            tmp = result.copy()
            str = word[start:pos + 1]
            print("entry", idx, start, pos, str)
            #new = pos - period * count + 1, period, count, idx
            new = pos - length * count + 1, length, count, idx
            # tmp.append((pos, idx, str))
            tmp.append(new)
            yield from walk2(inv, pmrs, max_cover, overlap, word, start - 1, tmp)

        # start = max(pmrs[idx][0], pos - pmrs[idx][1] * pmrs[idx][2] + 1)
        # value = pos - start + 1
        # print("value", value)
        # inc = max_cover[pos] - max_cover[start - 1]
        # print(f"max(pos): {max_cover[pos]} max(start - 1): {max_cover[start - 1]} inc: {inc}")
        #
        # tmp = result.copy()
        # if value + max_cover[start - 1] == max_cover[pos]:
        #     str = word[start:pos + 1]
        #     print("entry (whole)", idx, start, pos, str)
        #     tmp.append((pos, idx, str))
        #     walk2(inv, pmrs, max_cover, overlap, word, start - 1, tmp)
        #
        # tmp = result.copy()
        # if overlap[idx] and pos - overlap[idx] > 0:
        #     str = word[overlap[idx]:pos + 1]
        #     print("entry (no overlap conflict)", idx, overlap[idx], pos, str)
        #     tmp.append((pos, idx, str))
        #     walk2(inv, pmrs, max_cover, overlap, word, pos - 1, tmp)
        #
        # tmp = result.copy()
        # if overlap[idx] and pos - overlap[idx] == 0:
        #     print("No entry because of overlap")
        #     walk2(inv, pmrs, max_cover, overlap, word, pos - 1, tmp)
        #
        # tmp = result.copy()
        # if not overlap[idx]:
        #     start = pos - pmrs[idx][1] * pmrs[idx][2] + 1
        #     str = word[start:pos + 1]
        #     print("entry (no overlap)", idx, start, pos, str)
        #     tmp.append((pos, idx, str))
        #     walk2(inv, pmrs, max_cover, overlap, word, start - 1, tmp)



def fill(subs, word):
    desc = []

    end = -1
    total = 0
    for start, period, count, _ in subs:

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


def sol2path(solutions, max_cover):
    # print("sol2path")
    for start in range(len(max_cover) - 1, 0, -1):
        # print("s2p loop")
        if max_cover[start] != max_cover[-1]:
            # Don't look for paths without maximal result
            break
        yield from walk(solutions, max_cover[-1], start, [])


def sol2path3(inv, pmrs, max_cover, overlap, word):
    print("sol2path")
    for start in range(len(max_cover) - 1, 0, -1):
        print("s2p loop")
        if max_cover[start] != max_cover[-1]:
            # Don't look for paths without maximal result
            break
        yield from walk3(inv, pmrs, max_cover, overlap, word, start, [])


def path2hgvs(path, word):
    return ";".join(fill(path, word)[0])


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
    max_cover, solutions = cover(word, pmrs, hgvs=True)
    print_tables(n, word, inv, max_cover)
    # print(brute_cover(word, pmrs))

    print("solutions:")
    for idx, sol in enumerate(solutions):
        print(idx, sol)
    paths = sol2path(solutions, max_cover)
    print([path2hgvs(path, word) for path in paths])
    print("New")
    # paths = walk3(inv, pmrs, max_cover, overlap, word, len(word) - 1, [])
    # print(list(paths))
    # print([path2hgvs(path, word) for path in paths])
    paths = sol2path3(inv, pmrs, max_cover, overlap, word)
    print([path2hgvs(path, word) for path in paths])


if __name__ == "__main__":
    main()
