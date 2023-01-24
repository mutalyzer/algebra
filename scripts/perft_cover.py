from algebra.extractor.cover import cover_q, find_pmrs, inv_array
from algebra.extractor.words import non_isomorphic_binary_words
import math
import sys


def calculate_length(length):
    stat_work_min = math.inf
    stat_work_max = 0
    stat_work_sum = 0

    stat_pmrs_max = 0
    stat_pmrs_pos_max = 0
    stat_pmrs_sum = 0

    word_work_max = ""
    word_pmrs_max = ""
    word_pmrs_pos_max = ""

    count = 0
    for word in non_isomorphic_binary_words(length):
        count += 1
        pmrs = sorted(find_pmrs(word))
        inv = inv_array(length, pmrs)

        pmrs_counts = [len(x) for x in inv]
        if len(pmrs) >= stat_pmrs_max:
            stat_pmrs_max = len(pmrs)
            word_pmrs_max = word

        if max(pmrs_counts) >= stat_pmrs_pos_max:
            stat_pmrs_pos_max = max(pmrs_counts)
            word_pmrs_pos_max = word

        stat_pmrs_sum += sum(pmrs_counts)

        max_cover, work = cover_q(word, pmrs, inv)
        stat_work_min = min(stat_work_min, work)

        if work >= stat_work_max:
            stat_work_max = work
            word_work_max = word

        stat_work_sum += work

    print(length, stat_work_min, stat_work_max, stat_work_sum / count, stat_pmrs_max, stat_pmrs_pos_max, stat_pmrs_sum / count / length)
    print(word_work_max, word_pmrs_max, word_pmrs_pos_max)


def main():
    for length in range(1, int(sys.argv[1]) + 1):
        calculate_length(length)


if __name__ == "__main__":
    main()

