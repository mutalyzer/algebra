import itertools
from algebra.extractor.repeats import extract_repeats, subswithseqs
from algebra.extractor.cover import brute_cover, cover, inv_array


def main():
    counter = 0
    for length in range(2, 21):
        for word in ["".join(x) for x in itertools.product("AC", repeat=length)]:

            print(counter, word)

            pmrs = extract_repeats(word)

            if len(pmrs) < 2:
                continue

            print(pmrs)
            print(len(pmrs))

            inv = inv_array(len(word), pmrs)
            print(inv)
            if any([len(x) > 1 for x in inv]):
                continue

            subswithseqs(pmrs, word)

            cover_max = cover(word, pmrs)[-1]
            brute_max = brute_cover(word, pmrs)

            print(counter, word, cover_max, brute_max)
            assert cover_max == brute_max

            counter += 1
            print()


if __name__ == '__main__':
    main()
