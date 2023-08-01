import random
from algebra import compare
from algebra.relations.sequence_based import compare as compare_sequence
from algebra.relations.supremal_based import find_supremal, spanning_variant, compare as compare_supremal
from algebra.utils import random_sequence, random_variants
from algebra.variants import patch


def drive():
    reference = random_sequence(42, 1)
    lhs_var = list(random_variants(reference))
    lhs_seq = patch(reference, lhs_var)
    rhs_var = list(random_variants(reference))
    rhs_seq = patch(reference, rhs_var)

    lhs_span = spanning_variant(reference, lhs_seq, lhs_var)
    rhs_span = spanning_variant(reference, rhs_seq, rhs_var)
    lhs_sup, *_ = find_supremal(reference, lhs_span, offset=1)
    rhs_sup, *_ = find_supremal(reference, rhs_span, offset=1)

    rel_seq = compare_sequence(reference, lhs_seq, rhs_seq)
    rel_sup = compare_supremal(reference, lhs_sup, rhs_sup)
    rel_var = compare(reference, lhs_var, rhs_var)

    assert rel_seq == rel_sup == rel_var, (reference, lhs_var, rhs_var, rel_seq.value, rel_sup.value, rel_var.value)


def main():
    # random.seed(42)
    count = 0
    while True:
        drive()
        count += 1
        if count % 10_000 == 0:
            print(count)


if __name__ == "__main__":
    main()
