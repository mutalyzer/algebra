from argparse import ArgumentParser
from sys import stderr

from algebra.utils import fasta_sequence
from algebra.variants import parse_hgvs


def check_ins_ins(ncbi, alg):
    if (
        "ins" in ncbi
        and "del" not in ncbi
        and "ins" in alg
        and "del" not in alg
        and "[" in alg
    ):
        return True
    else:
        return False


def check_ins_repeat(ncbi, alg):
    if (
        "ins" in ncbi
        and "del" not in ncbi
        and "del" not in alg
        and "ins" not in alg
        and "[" in alg
    ):
        return True
    else:
        return False


def check_repeat_repeat(ncbi, alg, reference):
    if (
        "dup" not in ncbi
        and "del" not in ncbi
        and "ins" not in ncbi
        and "dup" not in alg
        and "del" not in alg
        and "ins" not in alg
        and ncbi.count("[") == 1
        and alg.count("[") == 1
        and parse_hgvs(ncbi, reference) == parse_hgvs(alg, reference)
    ):
        return True
    else:
        return False


def check_del_repeat(ncbi, alg):
    if (
        "del" in ncbi
        and "ins" not in ncbi
        and "del" not in alg
        and "ins" not in alg
        and "[" in alg
    ):
        return True
    else:
        return False


def check_dup_repeat(ncbi, alg):
    if "dup" in ncbi and "del" not in alg and "ins" not in alg and "[" in alg:
        return True
    else:
        return False


def main():
    parser = ArgumentParser(
        description="Compare the HGVS descriptions of algebra and NCBI"
    )
    parser.add_argument("--ncbi", required=True, type=str, help="Path to the NCBI file")
    parser.add_argument("--alg", required=True, type=str, help="Path to the algebra file")
    parser.add_argument("--fasta", required=True, type=str, help="Path to the fasta reference file")
    args = parser.parse_args()

    with open(args.fasta, encoding="utf-8") as file:
        reference = fasta_sequence(file)

    descriptions_ncbi = {}
    with open(args.ncbi, encoding="utf-8") as file:
        for line in file:
            description_spdi, description_hgvs = line.strip().split(" | ")
            descriptions_ncbi[description_spdi] = description_hgvs.split("g.")[1]

    descriptions_alg = {}
    with open(args.alg, encoding="utf-8") as file:
        for line in file:
            if "skipped" not in line:
                description_spdi, _, _, _, description_hgvs = line.strip().split(" ")
                descriptions_alg[description_spdi] = description_hgvs

    summary = {
        "common": {},
        "del_repeat": {},
        "dup_repeat": {},
        "ins_repeat": {},
        "ins_ins": {},
        "repeat_repeat": {},
        "other": {},
    }

    for description_spdi in descriptions_ncbi:
        if description_spdi in descriptions_alg:
            ncbi = descriptions_ncbi[description_spdi]
            alg = descriptions_alg[description_spdi]
            print(description_spdi, ncbi, alg)
            if ncbi == alg:
                summary["common"][description_spdi] = (ncbi, alg)
            else:
                print(description_spdi, ncbi, alg)
                if check_del_repeat(ncbi, alg):
                    summary["del_repeat"][description_spdi] = (ncbi, alg)
                elif check_dup_repeat(ncbi, alg):
                    summary["dup_repeat"][description_spdi] = (ncbi, alg)
                elif check_ins_repeat(ncbi, alg):
                    summary["ins_repeat"][description_spdi] = (ncbi, alg)
                elif check_ins_ins(ncbi, alg):
                    summary["ins_ins"][description_spdi] = (ncbi, alg)
                elif check_repeat_repeat(ncbi, alg, reference):
                    summary["repeat_repeat"][description_spdi] = (ncbi, alg)
                else:
                    summary["other"][description_spdi] = (ncbi, alg)

    for case in summary:
        print(f"{case}: {len(summary[case])}", file=stderr)
    for case in summary:
        if case != "common":
            print(f"\n{case}\n-----\n", file=stderr)
            for d in summary[case]:
                print(d, summary[case][d][0], summary[case][d][1], file=stderr)


if __name__ == "__main__":
    main()
