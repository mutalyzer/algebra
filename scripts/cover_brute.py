import itertools
from algebra.extractor.cover import cover, brute_cover_alt, cover_length, find_pmrs, overlapping, inv_array, cartesian_cover
from algebra.extractor.walker import path2hgvs, inv2paths, unique_pmrs


def main():
    counter = 0
    for length in range(21, 22):
        for word in ["".join(x) for x in itertools.product("AC", repeat=length)]:

            print(counter, word)

            pmrs = sorted(find_pmrs(word))

            if len(pmrs) < 2:
                continue

            # covers = list(brute_cover_alt(word, pmrs))
            # bmax = max(map(cover_length, covers))
            # brute_results = {path2hgvs(c, word) for c in covers if cover_length(c) == bmax}
            # print(brute_results)

            covers = list(cartesian_cover(pmrs))
            # print(word, covers)
            bmax = max(map(cover_length, covers))
            brute_results = {path2hgvs(c, word) for c in covers if cover_length(c) == bmax}
            print(brute_results)

            overlap = overlapping(pmrs)
            n = len(word)
            inv = inv_array(n, pmrs)

            max_cover = cover(word, pmrs)
            paths = inv2paths(inv, pmrs, max_cover, overlap, word)
            inv_results = {path2hgvs(path, word) for path in paths if unique_pmrs(path)}
            print(inv_results)

            assert brute_results == inv_results

            counter += 1


if __name__ == '__main__':
    main()
