import itertools
from algebra.extractor.repeats import brutepower, extract_repeats, subswithseqs
from algebra.extractor.cover import brute_cover, cover

def fib_word(n):
    fw0 = "A"
    if n <= 0:
        return fw0

    fw1 = "C"
    for _ in range(1, n):
        fw0, fw1 = fw1, fw1 + fw0
    return fw1


counter = 0
word = fib_word(9)

print(counter, word)

pmrs = extract_repeats(word)

print(pmrs)
print(len(pmrs))
subswithseqs(pmrs, word)

cover_max = cover(word, pmrs)[-1]
    #brute_max = brutepower(pmrs, word)[0]
    # brute_max = brute_cover(word, pmrs)
    #
    # print(counter, word, cover_max, brute_max)
    # assert cover_max == brute_max
    #
    # counter += 1

