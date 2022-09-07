TESTS = [
    ("AAB",
     [
        (0, 1, 2, 0),  # A
     ],
     [
        (1, 2),  # 0
     ],
     #  0  1  2
     #  A  A  B
     #  .  0  .
     #  0  2  2
     [[], [0], []],
     [None, 0, None],
     [0, 2, 2],
     "A[2];B",
    ),

    ("AABAB",
     [
        (0, 1, 2, 0),  # A
        (1, 2, 2, 0),  # AB
     ],
     [
        (1, 2),  # 0
        (4, 5),  # 1
     ],
     #  0  1  2  3  4
     #  A  A  B  A  B
     #  .  0  .  .  1
     #  0  2  2  2  4
     [[], [0], [], [], [1]],
     [None, 0, None, None, 1],
     [0, 2, 2, 2, 4],
     "A;AB[2]",
    ),

    ("AABC",
     [
        (0, 1, 2, 0),  # A
     ],
     [
        (1, 2),  # 0
     ],
     #  0  1  2  3
     #  A  A  B  C
     #  .  0  .  .
     #  0  2  2  2
     [[], [0], [], []],
     [None, 0, None, None],
     [0, 2, 2, 2],
     "A[2];BC",
    ),

    ("AABAAB",
     [
        (0, 1, 2, 0),  # A
        (0, 3, 2, 0),  # AAB
        (3, 1, 2, 0),  # A
     ],
     [
        (1, 2),  # 0
        (5, 6),  # 1
        (4, 5),  # 2
     ],
     #  0  1  2  3  4  5
     #  A  A  B  A  A  B
     #  .  0  .  .  2  1
     #  0  2  2  2  4  6
     [[], [0], [], [], [2], [1]],
     [None, 0, None, None, 2, 1],
     [0, 2, 2, 2, 4, 6],
     "AAB[2]",
    ),

    ("AABAABA",
     [
        (0, 1, 2, 0),  # A
        (0, 3, 2, 1),  # AAB
        (3, 1, 2, 0),  # A
     ],
     [
        (1, 2),  # 0
        (5, 7),  # 1
        (4, 5),  # 2
     ],
     #  0  1  2  3  4  5  6
     #  A  A  B  A  A  B  A
     #  .  0  .  .  2  1  1
     #  0  2  2  2  4  6  6
     [[], [0], [], [], [2], [1], [1]],
     [None, 0, None, None, 2, 1, 1],
     [0, 2, 2, 2, 4, 6, 6],
     "AAB[2];A",  # or "A;ABA[2]"
    ),

    ("AABAABB",
     [
        (0, 1, 2, 0),  # A
        (0, 3, 2, 0),  # AAB
        (3, 1, 2, 0),  # A
        (5, 1, 2, 0),  # B
     ],
     [
        (1, 2),  # 0
        (5, 6),  # 1
        (4, 5),  # 2
        (6, 7),  # 3
     ],
     #  0  1  2  3  4  5  6
     #  A  A  B  A  A  B  B
     #  .  0  .  .  2  1  3
     #  0  2  2  2  4  6  6
     [[], [0], [], [], [2], [1], [3]],
     [None, 0, None, None, 2, 1, 3],
     [0, 2, 2, 2, 4, 6, 6,],
     "AAB[2];B",  # or "A[2];B;A[2];B[2]"
    ),

    ("AABBAABB",
     [
        (0, 1, 2, 0),  # A
        (0, 4, 2, 0),  # AABB
        (2, 1, 2, 0),  # B
        (4, 1, 2, 0),  # A
        (6, 1, 2, 0),  # B
     ],
     [
        (1, 2),  # 0
        (7, 8),  # 1
        (3, 4),  # 2
        (5, 6),  # 3
        (7, 8),  # 4
     ],
     #  0  1  2  3  4  5  6  7
     #  A  A  B  B  A  A  B  B
     #  .  0  .  2  .  3  .  1
     #                       4
     #  0  2  2  4  4  6  6  8
     [[], [0], [], [2], [], [3], [], [1, 4]],
     [None, 0, None, 2, None, 3, None, 1],     # CONDITION I: IF only one, choose the larger unit
     [0, 2, 2, 4, 4, 6, 6, 8],
     "AABB[2]",  # or "A[2];B[2];A[2];B[2]"
    ),

    ("AAABAAAB",
     [
        (0, 1, 3, 0),  # A
        (0, 4, 2, 0),  # AAAB
        (4, 1, 3, 0),  # A
     ],
     [
        (1, 3),  # 0
        (7, 8),  # 1
        (5, 7),  # 2
     ],
     #  0  1  2  3  4  5  6  7
     #  A  A  A  B  A  A  A  B
     #  .  0  0  .  .  2  2  1
     #  0  2  3  3  3  5  6  8
     [[], [0], [0], [], [], [2], [2], [1]],
     [None, 0, 0, None, None, 2, 2, 1],
     [0, 2, 3, 3, 3, 5, 6, 8],
     "AAAB[2]",
    ),

    ("ACCACCTCT",
     [
        (0, 3, 2, 0),  # ACC
        (1, 1, 2, 0),  # C
        (4, 1, 2, 0),  # C
        (5, 2, 2, 0),  # CT
     ],
     [
        (5, 6),  # 0
        (2, 3),  # 1
        (5, 6),  # 2
        (8, 9),  # 3
     ],
     #  0  1  2  3  4  5  6  7  8
     #  A  C  C  A  C  C  T  C  T
     #  .  .  1  .  .  0  .  .  3
     #                 2
     #  0  0  2  2  2  6  6  6  6
     [[], [], [1], [], [], [0, 2], [], [], [3]],
     [None, None, 1, None, None, 0, None, None, 3],  # CONDITION I
     [0, 0, 2, 2, 2, 6, 6, 6, 6],
     "ACC[2];TCT",  # or "A;C[2];AC;CT[2]"
    ),

    ("TCCTCCACCA",
     [
        (0, 3, 2, 0),  # TCC
        (1, 1, 2, 0),  # C
        (4, 1, 2, 0),  # C
        (4, 3, 2, 0),  # CCA
        (7, 1, 2, 0),  # C
     ],
     [
        (5, 6),   # 0
        (2, 3),   # 1
        (5, 6),   # 2
        (9, 10),  # 3
        (8, 9),   # 4
     ],
     #  0  1  2  3  4  5  6  7  8  9
     #  T  C  C  T  C  C  A  C  C  A
     #  .  .  1  .  .  0  .  .  4  3
     #                 2
     #  0  0  2  2  2  6  6  6  8  8
     [[], [], [1], [], [], [0, 2], [], [], [4], [3]],
     [None, None, 1, None, None, 0, None, None, 4, 3],  # CONDITION I
     [0, 0, 2, 2, 2, 6, 6, 6, 8, 8],
     "TCC[2];A;C[2];A",  # or "T;C[2];T;CCA[2]"
    ),

    ("ABCABCABCDABCD",
     [
        (0, 3, 3, 0),  # ABC
        (6, 4, 2, 0),  # ABCD
     ],
     [
        (5, 9),    # 0
        (13, 14),  # 1
     ],
     #  0  1  2  3  4  5  6  7  8  9 10 11 12 13
     #  A  B  C  A  B  C  A  B  C  D  A  B  C  D
     #  .  .  .  .  .  0  0  0  0  .  .  .  .  1
     #  0  0  0  0  0  6  6  6  9  9  9  9  9 14
     [[], [], [], [], [], [0], [0], [0], [0], [], [], [], [], [1]],
     [None, None, None, None, None, 0, None, None, 0, None, None, None, None, 1],
     [0, 0, 0, 0, 0, 6, 6, 6, 9, 9, 9, 9, 9, 14],
     "ABC[2];ABCD[2]",
    ),

    ("CCACACAC",
     [
        (0, 1, 2, 0),  # C
        (1, 2, 3, 1),  # CA
     ],
     [
        (1, 2),  # 0
        (4, 8),  # 1
     ],
     #  0  1  2  3  4  5  6  7
     #  C  C  A  C  A  C  A  C
     #  .  0  .  .  1  1  1  1
     #  0  2  2  2  4  6  6  8
     [[], [0], [], [], [1], [1], [1], [1]],
     [None, 0, None, None, 1, 1, 1, 1],
     [0, 2, 2, 2, 4, 6, 6, 8],
     "C[2];AC[3]",
    ),

    ("CACATCACATCACATCACAT",
     [
        (0, 2, 2, 0),   # CA
        (0, 5, 4, 0),   # CACAT
        (5, 2, 2, 0),   # CA
        (10, 2, 2, 0),  # CA
        (15, 2, 2, 0),  # CA
     ],
     [
        (3, 4),    # 0
        (9, 20),   # 1
        (8, 9),    # 2
        (13, 14),  # 3
        (18, 19),  # 4
     ],
     #  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
     #  C  A  C  A  T  C  A  C  A  T  C  A  C  A  T  C  A  C  A  T
     #  .  .  .  0  .  .  .  .  2  1  1  1  1  1  1  1  1  1  1  1
     #                                         3              4
     #  0  0  0  4  4  4  4  4  8 10 10 10 10 14 15 15 15 15 19 20
     [[], [], [], [0], [], [], [], [], [2], [1], [1], [1], [1], [1, 3], [1], [1], [1], [1], [1, 4], [1]],
     [None, None, None, 0, None, None, None, None, 2, 1, None, None, None, 3, 1, None, None, None, 4, 1],
     [0, 0, 0, 4, 4, 4, 4, 4, 8, 10, 10, 10, 10, 14, 15, 15, 15, 15, 19, 20],
     "CACAT[4]"
    ),

    ("AACAACAAAACAACA",
     [
        (0, 1, 2, 0),   # A
        (0, 3, 2, 2),   # AAC
        (3, 1, 2, 0),   # A
        (3, 5, 2, 0),   # AACAA
        (6, 1, 4, 0),   # A
        (8, 3, 2, 1),   # AAC
        (11, 1, 2, 0),  # A
     ],
     [
        (1, 2),    # 0
        (5, 8),    # 1
        (4, 5),    # 2
        (12, 13),  # 3
        (7, 10),   # 4
        (13, 15),  # 5
        (12, 13),  # 6
     ],
     #  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
     #  A  A  C  A  A  C  A  A  A  A  C  A  A  C  A
     #  .  0  .  .  2  1  1  1  4  4  .  .  3  5  5
     #                       4              6
     #  0  2  2  2  4  6  6  8  9 10 10 10 12 14 15
     [[], [0], [], [], [2], [1], [1], [1, 4], [4], [4], [], [], [3, 6], [5], [5]],
     [None, 0, None, None, 2, 1, 1, 4, 4, 4, None, None, 3, 5, 5],   # CONDITION I + CONTIDITION II (overlap choose right)
     [0, 2, 2, 2, 4, 6, 6, 8, 9, 10, 10, 10, 12, 14, 15],
     "AAC[2];A[3];ACA[2]",
    ),

    ("AAAAAACAAAAAACAAAAAACACACAAAAAAA",
     [
        (0, 1, 6, 0),   # A
        (0, 7, 3, 1),   # AAAAAAC
        (7, 1, 6, 0),   # A
        (14, 1, 6, 0),  # A
        (19, 2, 3, 1),  # AC
        (25, 1, 7, 0),  # A
     ],
     [
        (1, 6),    # 0
        (13, 22),  # 1
        (8, 13),   # 2
        (15, 20),  # 3
        (22, 26),  # 4
        (26, 32),  # 5
     ],
     #  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
     #  A  A  A  A  A  A  C  A  A  A  A  A  A  C  A  A  A  A  A  A  C  A  C  A  C  A  A  A  A  A  A  A
     #  .  0  0  0  0  0  .  .  2  2  2  2  2  1  1  1  1  1  1  1  1  1  4  4  4  4  5  5  5  5  5  5
     #                                               3  3  3  3  3
     #  0  2  3  4  5  6  6  6  8  9 10 11 12 14 14 16 17 18 19 20 21 21 23 24 25 26 27 28 29 30 31 32
     [[], [0], [0], [0], [0], [0], [], [], [2], [2], [2], [2], [2], [1], [1], [1, 3], [1, 3], [1, 3], [1, 3], [1, 3], [1], [1], [4], [4], [4], [4], [5], [5], [5], [5], [5], [5]],
     [None, 0, 0, 0, 0, 0, None, None, 2, 2, 2, 2, 2, 1, 1, 3, 3, 3, 3, 3, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5],
     [0, 2, 3, 4, 5, 6, 6, 6, 8, 9, 10, 11, 12, 14, 14, 16, 17, 18, 19, 20, 21, 21, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32],
     "AAAAAAC[2];A[6];CA[3];A[6]",
    ),

    ("CATCATACATACTACTAAAAA",
     [
        (0, 3, 2, 0),   # CAT
        (3, 4, 2, 1),   # CATA
        (9, 3, 2, 2),   # TAC
        (16, 1, 5, 0),  # A
     ],
     [
        (5, 6),    # 0
        (10, 12),  # 1
        (14, 17),  # 2
        (17, 21),  # 3
     ],
     #  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
     #  C  A  T  C  A  T  A  C  A  T  A  C  T  A  C  T  A  A  A  A  A
     #  .  .  .  .  .  0  .  .  .  .  1  1  .  .  2  2  2  3  3  3  3
     #  0  0  0  0  0  6  6  6  6  6  8  8  8  8 12 12 14 14 16 17 18
     [[], [], [], [], [], [0], [], [], [], [], [1], [1], [], [], [2], [2], [2], [3], [3], [3], [3]],
     [None, None, None, None, None, 0, None, None, None, None, 1, 1, None, None, 2, 2, 2, 3, 3, 3, 3],
     [0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 8, 8, 8, 8, 12, 12, 14, 14, 16, 17, 18],
     "CAT;CATA[2];CTA[2];A[4]",
    ),

    ("CTACTACTATATA",
     [
        (0, 3, 3, 0),  # CTA
        (7, 2, 3, 0),  # TA
     ],
     [
         (5, 9),    # 0
         (10, 13),  # 1
     ],
     #  0  1  2  3  4  5  6  7  8  9 10 11 12
     #  C  T  A  C  T  A  C  T  A  T  A  T  A
     #  .  .  .  .  .  0  0  0  0  .  1  1  1
     #  0  0  0  0  0  6  6  6  9  9 10 10 13
     [[], [], [], [], [], [0], [0], [0], [0], [], [1], [1], [1]],
     [None, None, None, None, None, 0, None, None, 0, None, 1, None, 1],
     [0, 0, 0, 0, 0, 6, 6, 6, 9, 9, 10, 10, 13],
     "CTA[3];TA[2]",
    ),

    ("AATAAATAAAAA",
     [
        (0, 1, 2, 0),  # A
        (0, 4, 2, 2),  # AATA
        (3, 1, 3, 0),  # A
        (7, 1, 5, 0),  # A
     ],
     [
        (1, 2),   # 0
        (7, 10),  # 1
        (4, 6),   # 2
        (8, 12),  # 3
     ],
     #  0  1  2  3  4  5  6  7  8  9 10 11
     #  A  A  T  A  A  A  T  A  A  A  A  A
     #  .  0  .  .  2  2  .  1  1  1  3  3
     #                          3  3
     #  0  2  2  2  4  5  5  8  8 10 11 12
     [[], [0], [], [], [2], [2], [], [1], [1, 3], [1, 3], [3], [3]],
     [None, 0, None, None, 2, 2, None, 1, 3, 3, 3, 3],  # CONDITION II
     [0, 2, 2, 2, 4, 5, 5, 8, 8, 10, 11, 12],
     "AATA[2];A[4]",  # or "A[2];TAAA[2];A[2]"
    ),

    ("TAATAAAAAA",
     [
        (0, 3, 2, 0),  # TAA
        (1, 1, 2, 0),  # A
        (4, 1, 6, 0),  # A
     ],
     [
        (5, 6),   # 0
        (2, 3),   # 1
        (5, 10),  # 2
     ],
     #  0  1  2  3  4  5  6  7  8  9
     #  T  A  A  T  A  A  A  A  A  A
     #  .  .  1  .  .  0  2  2  2  2
     #                 2
     #  0  0  2  2  2  6  6  8  9 10
     [[], [], [1], [], [], [0, 2], [2], [2], [2], [2]],
     [None, None, 1, None, None, 0, 2, 2, 2, 2],  # CONDITION prefix choose left
     [0, 0, 2, 2, 2, 6, 6, 8, 9, 10],
     "TAA[2];A[4]",
    ),

    ("TACCACCCC",
     [
        (1, 3, 2, 0),  # ACC
        (2, 1, 2, 0),  # C
        (5, 1, 4, 0),  # C
     ],
     [
        (6, 7),  # 0
        (3, 4),  # 1
        (6, 9),  # 2
     ],
     #  0  1  2  3  4  5  6  7  8
     #  T  A  C  C  A  C  C  C  C
     #  .  .  .  1  .  .  0  2  2
     #                    2
     #  0  0  0  2  2  2  6  6  8
     [[], [], [], [1], [], [], [0, 2], [2], [2]],
     [None, None, None, 1, None, None, 0, 2, 2],  # CONDITION prefix
     [0, 0, 0, 2, 2, 2, 6, 6, 8],
     "T;ACC[2];C[2]",
    ),
]


def pmr_interval(pmr):
    start, period, count, remainder = pmr
    return start + 2 * period - 1, start + count * period + remainder


def inv_array(n, pmrs):
    inv = [[] for _ in range(n)]
    for idx, pmr in enumerate(pmrs):
        for pos in range(*pmr_interval(pmr)):
            inv[pos].append(idx)
    return inv


def inv_array_alt(n, pmrs):
    inv = [[] for _ in range(n)]
    for idx, pmr in enumerate(pmrs):
        start, period, count, remainder = pmr
        for c in range(1, count):
            for offset in range(remainder + 1):
                inv[start + (c + 1) * period + offset - 1].append(idx)
    return inv


def cover(word, pmrs, inv):
    n = len(word)
    max_cover = [0] * n

    hwm = -1
    for pos in range(1, n):
        # default to the previous value
        value = max_cover[pos - 1]

        if inv[pos] is not None:
            start, period, count, remainder = pmrs[inv[pos]]

            if hwm >= start:
                print("check hwm", hwm)
                real_count = (pos - hwm) // period
                if real_count > 1:
                    real_length = real_count * period
                    value = max(value, max_cover[hwm] + real_length)

            real_count = (pos - start + 1) // period
            if real_count > 1:
                real_length = real_count * period
                prev_value = 0
                if pos - real_length > -1:
                    prev_value = max_cover[pos - real_length]

                if prev_value + real_length > value:
                    value = prev_value + real_length
                    if pos == start + period * count + remainder - 1:
                        print("update hwm", pos, pmrs[inv[pos]])
                        hwm = pos

        max_cover[pos] = value

    return max_cover


def main():
    for word, pmrs, pmr_ints, inv, inv_flat, max_cover, hgvs in TESTS:
        print(word)
        for idx, pmr in enumerate(pmrs):
            assert pmr_ints[idx] == pmr_interval(pmr)
        assert inv == inv_array(len(word), pmrs)

        #print(inv_array_alt(len(word), pmrs))

        test = cover(word, pmrs, inv_flat)
        if max_cover != test:
            print(test)
            print(max_cover)


if __name__ == "__main__":
    main()
