def cover(word, runs, n=0, repeats=[]):
    if n >= len(runs):
        def hgvs():
            start = 0
            for repeat in repeats:
                if repeat[2] == 1:
                    continue
                if repeat[0] > start:
                    yield word[start:repeat[0]]
                yield f"{word[repeat[0]:repeat[0] + repeat[1]]}[{repeat[2]}]"
                start = repeat[0] + repeat[1] * repeat[2]
            if start < len(word):
                yield word[start:]
        yield ";".join(hgvs())

    for idx, run in enumerate(runs[n:]):
        prev = (-1, 0, 0, 0) if not repeats else repeats[-1]
        for rotation in range(run[3] + 1):
            curr = run[0] + rotation, run[1], run[2], run[3] - rotation
            if prev[0] + prev[1] * prev[2] > curr[0]:
                curr = run[0] + rotation + run[1], run[1], run[2] - 1, run[3] - rotation

            repeats.append(curr)
            yield from cover(word, runs, n + idx + 1, repeats)
            repeats.pop()


def main():
    word = "CATCATACATACTACTAAAAA"
    runs = [
        ( 0, 3, 2, 0),
        ( 3, 4, 2, 1),
        ( 9, 3, 2, 2),
        (16, 1, 5, 0),
    ]
    for hgvs in set(cover(word, runs)):
        print(hgvs)


if __name__ == "__main__":
    main()
