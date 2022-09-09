from itertools import chain, combinations
from os.path import commonprefix
import sys
from algebra.utils import random_sequence
from algebra.extractor.repeats import brutepower, extract_repeats, subswithseqs
from algebra.extractor.cover import cover


def main():
    counter = 0
    while True:
        word = random_sequence(10, 5)
        pmrs = extract_repeats(word)
        if not pmrs:
            continue

        print(word)
        print(pmrs)
        print(len(pmrs))
        subswithseqs(pmrs, word)

        cover_max = cover(word, pmrs)[-1]
        brute_max = brutepower(pmrs, word)[0]

        print(counter, word, cover_max, brute_max)
        assert cover_max == brute_max
        counter += 1


if __name__ == '__main__':
    main()
