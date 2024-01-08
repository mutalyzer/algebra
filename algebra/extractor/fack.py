import sys


def fac(k):
    if k < 2:
        return

    yield (k,)
    if k < 4:
        return

    for i in range(2, k - 1):
        low = list(fac(i))
        high = list(fac(k - i))
        for left in low:
            for right in high:
                yield (*left, *right)


def main():
    for k in range(2, int(sys.argv[1])):
        factors = sorted(set(fac(k)))
        print(k, len(factors), factors)


if __name__ == "__main__":
    main()

