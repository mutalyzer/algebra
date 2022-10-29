from itertools import islice, product


def non_isomorphic_binary_words(length, alphabet="AC"):
    """Generates all non-isomorphic words of a certain length over a
    binary alphabet."""
    if length <= 0:
        yield ""
        return
    for suffix in product(alphabet[:2], repeat=length-1):
        yield alphabet[0] + "".join(suffix)


def non_isomorphic_binary_words_upto(length, alphabet="AC"):
    """Generates all non-isomorphic words up to a certain length over a
    binary alphabet."""
    for i in range(length + 1):
        yield from non_isomorphic_binary_words(i, alphabet)


def nth_fibonacci_word(n, alphabet="AC"):
    """Returns the n-th (finite) Fibonacci word."""
    if n < 0:
        return ""
    word0, word1 = alphabet[0], alphabet[1]
    for _ in range(n):
        word0, word1 = word0 + word1, word0
    return word0


def fibonacci_word(alphabet="AC"):
    """Generates the infinite Fibonacci word."""
    idx = 0
    n = 0
    word = nth_fibonacci_word(n)
    while True:
        try:
            yield word[idx]
            idx += 1
        except IndexError:
            n += 1
            word = nth_fibonacci_word(n)


def fibonacci_word_upto(length, alphabet="AC"):
    """Returns the prefix of a certain length of the infinite Fibonacci
    word."""
    if length < 0:
        return ""
    return "".join(islice(fibonacci_word(), length))


def lyndon_words(length, alphabet="AC"):
    """Generates all Lyndon words of a certain length over a ordered
    alphabet using Duval's algorithm."""
    w = [-1]
    while w:
        w[-1] += 1
        m = len(w)
        if m == length:
            yield "".join(alphabet[i] for i in w)

        while len(w) < length:
            w.append(w[-m])

        while w and w[-1] == len(alphabet) - 1:
            w.pop()
