from itertools import chain, combinations
from os.path import commonprefix
import sys
from algebra.utils import random_sequence


def main():
    while True:
        word = random_sequence(10, 5)
        pmrs = extract_repeats(word)
        print(pmrs)
        print(len(pmrs))
        subswithseqs(pmrs, word)

        koffer_max = koffer(pmrs, word)
        brute_max = brutepower(pmrs, word)

        print(word, koffer_max, brute_max)
        assert koffer_max == brute_max


if __name__ == '__main__':
    main()
