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

    args = parser.parse_args()
    if args.sequences:
        normal = list(extract(*args.sequences))
        dfa = list(extract_dfa(*args.sequences))
        if normal == dfa:
            print(f"equal: {[e.to_hgvs() for e in normal]}")
        else:
            print(f"- normal:\n{normal}")
            print(f"- dfa:\n{dfa}")

    else:
        count = 0
        while True:
            count += 1
            ref = random_sequence(100, 30)
            obs = random_sequence(100, 30)
            print(f"\n===\n{ref} {obs}\n---\n{count}")
            assert list(extract(ref, obs)) == list(extract_dfa(ref, obs))


if __name__ == "__main__":
    main()
