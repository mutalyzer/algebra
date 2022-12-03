import itertools
from algebra.extractor.cover import cover, brute_cover_alt, cover_length, find_pmrs, overlapping, inv_array, cartesian_cover
from algebra.extractor.walker import path2hgvs, unique_pmrs, walk


def main():
    counter = 0
    for length in range(2, 21):
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
            brute_results = sorted([path2hgvs(c, word) for c in covers if cover_length(c) == bmax])
            print(brute_results)

            overlap = overlapping(pmrs)
            n = len(word)
            inv = inv_array(n, pmrs)

            max_cover = cover(word, pmrs)
            paths = list(walk(inv, pmrs, max_cover, overlap, len(max_cover) - 1, []))
            all_max = max(map(cover_length, paths))
            paths_uniq = [path for path in paths if unique_pmrs(path)]
            uniq_max = max(map(cover_length, paths_uniq))
            inv_results_uniq = sorted([path2hgvs(path, word) for path in paths if unique_pmrs(path) and cover_length(path) == uniq_max])
            print(inv_results_uniq)

            assert brute_results == inv_results_uniq
            assert bmax == all_max == uniq_max == max_cover[-1]

            counter += 1


if __name__ == '__main__':
    main()
