from argparse import ArgumentParser

from algebra.extractor.extractor import extract
from algebra.extractor.extractor_dfa import extract as extract_dfa
from algebra.utils import random_sequence


def main():
    parser = ArgumentParser(
        "Check and compare the extracted description with the DFA approach"
    )
    parser.add_argument(
        "-s",
        "--sequences",
        action="store",
        nargs=2,
        metavar=("ref", "obs"),
        help="check for the specified sequences",
    )
    parser.add_argument(
        "-i",
        "--interval",
        action="store",
        nargs=2,
        metavar=("max", "min"),
        type=int,
        help="sequences length",
    )

    args = parser.parse_args()
    if args.sequences:
        normal = list(extract(*args.sequences))
        dfa = list(extract_dfa(*args.sequences))
        if normal == dfa:
            print(f"equal: {[e.to_hgvs() for e in normal]}")
        else:
            print(f"- normal:\n{[e.to_hgvs() for e in normal]}")
            print(f"- dfa:\n{[e.to_hgvs() for e in dfa]}")

    else:
        count = 0
        len_max = 10 if not args.interval else max(args.interval)
        len_min = 0 if not args.interval else min(args.interval)
        while True:
            count += 1
            ref = random_sequence(len_max, len_min)
            obs = random_sequence(len_max, len_min)
            print(f"\n=== [{len_max}, {len_min}]\n{ref} {obs}\n{len(ref)} {len(obs)}\n---\n{count}")
            assert list(extract(ref, obs)) == list(extract_dfa(ref, obs))


if __name__ == "__main__":
    main()
