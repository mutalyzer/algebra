from algebra.relations.sequence_based import compare as compare_sequence
from algebra.relations.supremal_based import find_supremal, spanning_variant, compare as compare_supremal
from algebra.utils import random_sequence, random_variants
from algebra.variants import patch
import random


def drive():
    reference = random_sequence(42, 1)
    lhs_var = list(random_variants(reference))
    lhs_seq = patch(reference, lhs_var)
    rhs_var = list(random_variants(reference))
    rhs_seq = patch(reference, rhs_var)

    try:
        lhs_span = spanning_variant(reference, lhs_seq, lhs_var)
        rhs_span = spanning_variant(reference, rhs_seq, rhs_var)
        lhs_sup = find_supremal(reference, lhs_span, offset=1)
        rhs_sup = find_supremal(reference, rhs_span, offset=1)
    except ValueError:
        return

    rel_seq = compare_sequence(reference, lhs_seq, rhs_seq)
    rel_sup = compare_supremal(reference, lhs_sup, rhs_sup)

    assert rel_seq == rel_sup, (reference, lhs_var, rhs_var, rel_seq.value, rel_sup.value)


def main():
    count = 0
    while True:
        drive()
        count += 1
        if count % 100000 == 0:
            print(count)


if __name__ == '__main__':
    # random.seed(42)
    main()
