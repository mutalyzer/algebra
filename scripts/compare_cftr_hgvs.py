import argparse


def main():
    parser = argparse.ArgumentParser(description="Compare the HGVS descriptions of algebra and NCBI")
    parser.add_argument("--refseq-id", required=True, type=str, help="RefSeq identifier")
    parser.add_argument("--ncbi", required=True, type=str, help="Path to the NCBI file")
    parser.add_argument("--alg", required=True, type=str, help="Path to the algebra file")

    args = parser.parse_args()

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

    # print(len(descriptions_ncbi))
    # print(len(descriptions_alg))
    common = 0
    del_repeat = 0
    dup_repeat = 0
    other = 0
    for description_spdi in descriptions_ncbi:
        if description_spdi in descriptions_alg:
            ncbi = descriptions_ncbi[description_spdi]
            alg = descriptions_alg[description_spdi]
            if ncbi == alg:
                common += 1
            elif "del" in ncbi and "ins" not in ncbi and "del" not in alg and "ins" not in alg and "[" in alg:
                # print(description_spdi, ncbi, alg)
                del_repeat += 1
            elif "dup" in ncbi and "del" not in alg and "ins" not in alg and "[" in alg:
                # print(description_spdi, ncbi, alg)
                dup_repeat += 1
            else:
                print(description_spdi, ncbi, alg)
                other += 1
    print("\n--------\n")
    print(f"common: {common}")
    print(f"del_repeat: {del_repeat}")
    print(f"dup_repeat: {dup_repeat}")
    print(f"other: {other}")


if __name__ == "__main__":
    main()
