from itertools import product
from multiprocessing import Pool
import sys

from algebra.extractor.cover import brute_cover, cover, find_pmrs
from algebra.extractor.perf import tic, toc
from algebra.extractor.words import non_isomorphic_binary_words_upto


def check_with_brute(word):
    pmrs = sorted(find_pmrs(word))
    max_cover = cover(word, pmrs)
    max_brute_cover = brute_cover(word, pmrs)
    if max_cover and max_brute_cover != max_cover[-1]:
        print(word, max_cover, max_brute_cover)
        assert False


def main():
    if len(sys.argv) < 2:
        print(f"usage: {sys.argv[0]} length", file=sys.stderr)
        return

    length = int(sys.argv[1])

    tic()
    with Pool() as pool:
        pool.map(check_with_brute, non_isomorphic_binary_words_upto(length))
    #list(map(check_with_brute, non_isomorphic_binary_words_upto(length)))
    toc()


if __name__ == "__main__":
    main()
