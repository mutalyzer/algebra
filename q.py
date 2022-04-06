import sys
from algebra.lcs.onp import edit as edit_onp
from algebra.lcs.wupp import edit as edit_test, lcs_graph
from algebra.lcs.efficient import edit as edit_gold
from pprint import pprint
import random


def compare(test, gold, f):
    assert len(test) == len(gold)
    assert len(test[0]) == len(gold[0])
    row_len = len(test)
    col_len = len(test[0])

    for row in range(row_len):
        for col in range(col_len):
            if gold[row][col] is not None and gold[row][col] <= f:
                assert gold[row][col] == test[row][col]


def main():
    min_rand = 1
    max_rand = 42

    if len(sys.argv) == 1:
        reference = "".join(random.choice("ACGT") for _ in range(random.randint(min_rand, max_rand)))
        observed = "".join(random.choice("ACGT") for _ in range(random.randint(min_rand, max_rand)))
    elif len(sys.argv) == 3:
        reference = sys.argv[1]
        observed = sys.argv[2]
    else:
        raise("usage")

    print(reference, observed)

    dist_test, matrix_test, lcs_nodes = edit_test(reference, observed)
    dist_gold, _, matrix_gold = edit_gold(reference, observed)
    print(dist_test)
    # pprint(matrix_test)
    print(dist_gold)
    # pprint(matrix_gold)

    assert dist_test == dist_gold

    compare(matrix_test, matrix_gold, dist_test)

    for level in lcs_nodes:
        print(level)

    lcs_graph(reference, observed, lcs_nodes)


if __name__ == '__main__':
    main()
