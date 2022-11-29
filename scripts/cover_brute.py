import itertools
from algebra.extractor.repeats import brutepower, subswithseqs
from algebra.extractor.cover import cover, brute_cover_alt, cover_length, find_pmrs, overlapping, inv_array, cartesian_cover
from algebra.extractor.walker import path2hgvs, inv2paths, sol2paths


def unique_pmrs(cover):
    # print(cover)
    pmrs = [x[3] for x in cover]
    # print(pmrs)
    if sorted(pmrs) == sorted(set(pmrs)):
        return True
    return False


def main():
    counter = 0
    for length in range(20, 30):
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
            # brute_results = {path2hgvs(c, word) for c in covers}
            # brute_results = {path2hgvs(c, word) for c in covers if cover_length(c) == bmax}
            brute_results = {path2hgvs(c, word) for c in covers if cover_length(c) == bmax and unique_pmrs(c)}
            print(brute_results)

            overlap = overlapping(pmrs)
            n = len(word)
            inv = inv_array(n, pmrs)

            max_cover, solutions = cover(word, pmrs, hgvs=True)

            # paths = sol2paths(word, solutions, max_cover)
            # sol_results = {path2hgvs(path, word) for path in paths}
            # print(sol_results)

            paths = inv2paths(inv, pmrs, max_cover, overlap, word)
            inv_results = {path2hgvs(path, word) for path in paths if unique_pmrs(path)}
            print(inv_results)

            # assert brute_results == sol_results == inv_results
            # assert sol_results == inv_results
            # assert brute_results == sol_results
            assert brute_results == inv_results

            counter += 1


if __name__ == '__main__':
    main()
