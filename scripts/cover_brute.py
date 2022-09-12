from itertools import chain, combinations
from os.path import commonprefix
import subprocess
import sys
from algebra.utils import random_sequence
from algebra.extractor.repeats import brutepower, extract_repeats, subswithseqs
from algebra.extractor.cover import cover
from algebra.extractor.repeats_alt import _repeats


def main():
    counter = 0
    while True:
        word = random_sequence(10, 5)
        print(word)
        # pmrs1 = extract_repeats(word)
        pmrs = sorted([(x[1], len(x[0]), x[2], x[3]) for x in _repeats(word)])
        # print(pmrs1, pmrs)

        #p = subprocess.run(["../repeats/mreps/mreps", "-allowsmall", "-s", word], capture_output=True)

        #pmrs_mreps = eval(p.stdout)

        #assert pmrs == pmrs_mreps

        if not pmrs:
            continue

        # continue

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
