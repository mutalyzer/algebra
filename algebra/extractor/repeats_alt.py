from argparse import ArgumentParser


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


if __name__ == "__main__":
    main()
