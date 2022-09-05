from argparse import ArgumentParser

import pytest


def _window_fixed_size(seq, window_sz):
    repeats = []
    window_start = 0
    window_seq = seq[window_start : window_start + window_sz]
    idx = 0
    counts = 1
    while window_start + counts * window_sz + idx < len(seq):
        if window_seq[idx] == seq[window_start + counts * window_sz + idx]:
            idx += 1
            if idx == window_sz:
                counts += 1
                idx = 0
        else:
            if counts > 1:
                repeats.append((window_seq, window_start, counts, idx))
            window_start += 1
            window_seq = seq[window_start : window_start + window_sz]
            counts = 1
            idx = 0
    if counts > 1:
        repeats.append((window_seq, window_start, counts, idx))

    return repeats


def _repeats(seq):
    def _already_in(new_interval):
        for in_interval in intervals:
            if new_interval[0] >= in_interval[0] and new_interval[1] <= in_interval[1]:
                return True
        return False

    repeats = []
    intervals = []

    for window_size in range(1, len(seq) // 2 + 1):
        for repeat in _window_fixed_size(seq, window_size):
            interval = repeat[1], repeat[1] + len(repeat[0]) * repeat[2] + repeat[3]
            if _already_in(interval):
                continue
            else:
                intervals.append(interval)
                repeats.append(repeat)
    return repeats


def main():
    parser = ArgumentParser()
    parser.add_argument("seq", type=str, help="Sequence")
    args = parser.parse_args()

    for repeat in _repeats(args.seq):
        print(repeat)


@pytest.mark.parametrize(
    "sequence, repeats",
    [
        (
            "CATCATACATACTACTAAAAA",
            [
                ("A", 16, 5, 0),
                ("CAT", 0, 2, 0),
                ("TAC", 9, 2, 2),
                ("CATA", 3, 2, 1),
            ],
        ),
        (
            "ACCACCTCT",
            [
                ("C", 1, 2, 0),
                ("C", 4, 2, 0),
                ("CT", 5, 2, 0),
                ("ACC", 0, 2, 0),
            ],
        ),
        (
            "CCACACAC",
            [
                ("C", 0, 2, 0),
                ("CA", 1, 3, 1),
            ],
        ),
        (
            "CACATCACATCACATCACAT",
            [
                ("CA", 0, 2, 0),
                ("CA", 5, 2, 0),
                ("CA", 10, 2, 0),
                ("CA", 15, 2, 0),
                ("CACAT", 0, 4, 0),
            ],
        ),
        (
            "AACAACAAAACAACA",
            [
                ("A", 0, 2, 0),
                ("A", 3, 2, 0),
                ("A", 6, 4, 0),
                ("A", 11, 2, 0),
                ("AAC", 0, 2, 2),
                ("AAC", 8, 2, 1),
                ("AACAA", 3, 2, 0),
            ],
        ),
        (
            "AAAAAACAAAAAACAAAAAACACACAAAAAAA",
            [
                ("A", 0, 6, 0),
                ("A", 7, 6, 0),
                ("A", 14, 6, 0),
                ("A", 25, 7, 0),
                ("AC", 19, 3, 1),
                ("AAAAAAC", 0, 3, 1),
            ],
        ),
    ],
)
def test_repeats(sequence, repeats):
    assert repeats == _repeats(sequence)


if __name__ == "__main__":
    main()
