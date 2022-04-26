import random
import sys
from algebra.relations_efficient import compare as compare_gold
from algebra.relations import compare as compare_test


def main():
    min_rand = 0
    max_rand = 42

    if len(sys.argv) == 1:
        reference = "".join(random.choice("ACGT") for _ in range(random.randint(min_rand, max_rand)))
        lhs = "".join(random.choice("ACGT") for _ in range(random.randint(min_rand, max_rand)))
        rhs = "".join(random.choice("ACGT") for _ in range(random.randint(min_rand, max_rand)))
    elif len(sys.argv) == 4:
        reference = sys.argv[1]
        lhs = sys.argv[2]
        rhs = sys.argv[3]
    else:
        raise("usage")

    print(reference, lhs, rhs)

    assert compare_gold(reference, lhs, rhs) == compare_test(reference, lhs, rhs)


if __name__ == '__main__':
    main()
