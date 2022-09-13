import os
import subprocess
from algebra.utils import random_sequence
from algebra.extractor.repeats import extract_repeats, subswithseqs
from algebra.extractor.repeats_alt import _repeats


def main():
    counter = 0
    while True:
        word = random_sequence(42, 5)
        print(counter, word)

        # "Original"
        pmrs = extract_repeats(word)

        # "Mihai"
        pmrs_alt = sorted([(x[1], len(x[0]), x[2], x[3]) for x in _repeats(word)])

        # mreps if available
        mreps_bin = "../repeats/mreps/mreps"
        if os.path.isfile(mreps_bin) and os.access(mreps_bin, os.X_OK):
            p = subprocess.run([mreps_bin, "-allowsmall", "-s", word], capture_output=True)
            pmrs_mreps = eval(p.stdout)
            assert pmrs == pmrs_alt == pmrs_mreps
        else:
            assert pmrs == pmrs_alt

        if not pmrs:
            continue

        counter += 1

        #print(pmrs)
        #print(len(pmrs))
        #subswithseqs(pmrs, word)


if __name__ == '__main__':
    main()
