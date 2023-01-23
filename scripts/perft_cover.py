from algebra.extractor.cover import cover_q, find_pmrs
from algebra.extractor.words import non_isomorphic_binary_words
import math
import sys


def calculate_length(length):
    stat_work_min = math.inf
    stat_work_max = 0
    stat_work_sum = 0
    stat_work_count = 0

    for word in non_isomorphic_binary_words(length):
        pmrs = sorted(find_pmrs(word))
        max_cover, work = cover_q(word, pmrs)
        stat_work_min = min(stat_work_min, work)
        stat_work_max = max(stat_work_max, work)
        stat_work_count += 1
        stat_work_sum += work

    print(length, stat_work_min, stat_work_max, stat_work_sum / stat_work_count)


def main():
    for length in range(1, int(sys.argv[1]) + 1):
        calculate_length(length)


if __name__ == "__main__":
    main()

