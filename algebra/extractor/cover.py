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


def cover(word, runs, n=0, repeats=[]):
    if n >= len(runs):
        yield repeats

    for idx, run in enumerate(runs[n:]):
        prev = (-1, 0, 0, 0) if not repeats else repeats[-1]
        for rotation in range(run[3] + 1):
            curr = run[0] + rotation, run[1], run[2], run[3] - rotation
            if prev[0] + prev[1] * prev[2] > curr[0]:
                curr = run[0] + rotation + run[1], run[1], run[2] - 1, run[3] - rotation
                if curr[2] == 1:
                    continue

            repeats.append(curr)
            yield from cover(word, runs, n + idx + 1, repeats)
            repeats.pop()


def cover2len(repeats):
    return sum([x[1] * x[2] for x in repeats])


def maxcover(word, runs):
    return max([cover2len(x) for x in cover(word, runs)])


def koffer(word, runs):

    a = [None] * len(word)
    for idx, run in enumerate(runs):
        start = run[0] + run[1] * 2 - 1
        end = run[0] + run[1] * run[2] + run[3]
        print(run)
        print(start, end)

        for i in range(start, end):
            a[i] = idx

    print("run indices:", a)

    b = [0] * len(word)

    for idx, run_idx in enumerate(a):
        print("idx", idx)

        prev = 0
        if idx > 0:
            prev = b[idx - 1]

        q = 0
        if run_idx is not None:
            run = runs[run_idx]

            l = run[1] * ((idx - run[0] + 1) // run[1])
            print("l", l)

            k = (idx - run[0] + 1) % run[1]
            print("k", k)

            if idx > run[0] + 1:
                print("idx > run")
                q = l + b[run[0] + k]
                print(q)


        b[idx] = max(prev, q)
        print(f"new b[{idx}]: {b[idx]}")


    print(b)
    return a[-1]


def main():
    word = "CATCATACATACTACTAAAAA"
    # word = "BAABAA"
    runs = [
        # (0, 3, 2, 0),
        # (1, 1, 2, 0),
        # (4, 1, 2, 0),
        ( 0, 3, 2, 0),
        ( 3, 4, 2, 1),
        ( 9, 3, 2, 2),
        (16, 1, 5, 0),
    ]
    # for repeats in cover(word, runs):
    #     print(to_hgvs(word, repeats))
    #     print(cover2len(repeats))

    print(maxcover(word, runs))

    koffer(word, runs)


if __name__ == "__main__":
    main()
