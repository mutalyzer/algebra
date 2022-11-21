from itertools import product
import sys

from algebra.extractor.cover import brute_cover, cover, find_pmrs
from algebra.extractor.words import non_isomorphic_binary_words


def check_with_brute(word):
    pmrs = sorted(find_pmrs(word))
    max_cover = cover(word, pmrs)
    max_brute_cover = brute_cover(word, pmrs)
    if max_cover and max_brute_cover != max_cover[-1]:
        print(word, max_cover, max_brute_cover)
        assert False
    print(len(word), word, len(pmrs), max_brute_cover)


def main():
    if len(sys.argv) < 2:
        print(f"usage: {sys.argv[0]} length", file=sys.stderr)
        return

    length = int(sys.argv[1])

    for word in non_isomorphic_binary_words(length):
        check_with_brute(word)


if __name__ == "__main__":
    main()
