import pytest
from algebra.extractor.cover import brute_cover, cover, find_pmrs, inv_array, overlapping
from algebra.extractor.repeats import extract_repeats
from algebra.extractor.repeats_alt import _repeats as repeats_alt
from algebra.extractor.walker import path2hgvs, sol2path, sol2path3


TESTS = [
    ("AAB",
     [
         (0, 1, 2, 0),  # A
     ],
     #  0  1  2
     #  A  A  B
     #  .  0  .
     #  0  2  2
     [[], [0], []],  # flat
     [0, 2, 2],
     {"A[2];B"},
     ),

    ("AABAB",
     [
         (0, 1, 2, 0),  # A
         (1, 2, 2, 0),  # AB
     ],
     #  0  1  2  3  4
     #  A  A  B  A  B
     #  .  0  .  .  1
     #  0  2  2  2  4
     [[], [0], [], [], [1]],  # flat
     [0, 2, 2, 2, 4],
     {"A;AB[2]"},
     ),

    ("AABC",
     [
         (0, 1, 2, 0),  # A
     ],
     #  0  1  2  3
     #  A  A  B  C
     #  .  0  .  .
     #  0  2  2  2
     [[], [0], [], []],  # flat
     [0, 2, 2, 2],
     {"A[2];BC"},
     ),

    ("AABAAB",
     [
         (0, 1, 2, 0),  # A
         (0, 3, 2, 0),  # AAB
         (3, 1, 2, 0),  # A
     ],
     #  0  1  2  3  4  5
     #  A  A  B  A  A  B
     #  .  0  .  .  2  1
     #  0  2  2  2  4  6
     [[], [0], [], [], [2], [1]],  # flat
     [0, 2, 2, 2, 4, 6],
     {"AAB[2]"},
     ),

    ("AABAABA",
     [
         (0, 1, 2, 0),  # A
         (0, 3, 2, 1),  # AAB
         (3, 1, 2, 0),  # A
     ],
     #  0  1  2  3  4  5  6
     #  A  A  B  A  A  B  A
     #  .  0  .  .  2  1  1
     #  0  2  2  2  4  6  6
     [[], [0], [], [], [2], [1], [1]],  # flat
     [0, 2, 2, 2, 4, 6, 6],
     {
         "AAB[2];A",
         "A;ABA[2]",
     }
     ),

    ("AABAABB",
     [
         (0, 1, 2, 0),  # A
         (0, 3, 2, 0),  # AAB
         (3, 1, 2, 0),  # A
         (5, 1, 2, 0),  # B
     ],
     #  0  1  2  3  4  5  6
     #  A  A  B  A  A  B  B
     #  .  0  .  .  2  1  3
     #  0  2  2  2  4  6  6
     [[], [0], [], [], [2], [1], [3]],  # flat
     [0, 2, 2, 2, 4, 6, 6],
     {
         "AAB[2];B",
         "A[2];B;A[2];B[2]",
     }
     ),

    ("AABBAABB",
     [
         (0, 1, 2, 0),  # A
         (0, 4, 2, 0),  # AABB
         (2, 1, 2, 0),  # B
         (4, 1, 2, 0),  # A
         (6, 1, 2, 0),  # B
     ],
     #  0  1  2  3  4  5  6  7
     #  A  A  B  B  A  A  B  B
     #  .  0  .  2  .  3  .  1
     #                       4
     #  0  2  2  4  4  6  6  8
     [[], [0], [], [2], [], [3], [], [1, 4]],  # prefix
     [0, 2, 2, 4, 4, 6, 6, 8],
     {
        "AABB[2]",
        "A[2];B[2];A[2];B[2]",
     }
     ),

    ("AAABAAAB",
     [
         (0, 1, 3, 0),  # A
         (0, 4, 2, 0),  # AAAB
         (4, 1, 3, 0),  # A
     ],
     #  0  1  2  3  4  5  6  7
     #  A  A  A  B  A  A  A  B
     #  .  0  0  .  .  2  2  1
     #  0  2  3  3  3  5  6  8
     [[], [0], [0], [], [], [2], [2], [1]],  # flat
     [0, 2, 3, 3, 3, 5, 6, 8],
     {"AAAB[2]"},
     ),

    ("ACCACCTCT",
     [
         (0, 3, 2, 0),  # ACC
         (1, 1, 2, 0),  # C
         (4, 1, 2, 0),  # C
         (5, 2, 2, 0),  # CT
     ],
     #  0  1  2  3  4  5  6  7  8
     #  A  C  C  A  C  C  T  C  T
     #  .  .  1  .  .  0  .  .  3
     #                 2
     #  0  0  2  2  2  6  6  6  6
     [[], [], [1], [], [], [0, 2], [], [], [3]],  # prefix
     [0, 0, 2, 2, 2, 6, 6, 6, 6],
     {
         "ACC[2];TCT",
         "A;C[2];AC;CT[2]",
     }
     ),

    ("TCCTCCACCA",
     [
         (0, 3, 2, 0),  # TCC
         (1, 1, 2, 0),  # C
         (4, 1, 2, 0),  # C
         (4, 3, 2, 0),  # CCA
         (7, 1, 2, 0),  # C
     ],
     #  0  1  2  3  4  5  6  7  8  9
     #  T  C  C  T  C  C  A  C  C  A
     #  .  .  1  .  .  0  .  .  4  3
     #                 2
     #  0  0  2  2  2  6  6  6  8  8
     [[], [], [1], [], [], [0, 2], [], [], [4], [3]],  # prefix
     [0, 0, 2, 2, 2, 6, 6, 6, 8, 8],
     {
         "TCC[2];A;C[2];A",
         "T;C[2];T;CCA[2]",
     }
     ),

    ("ABCABCABCDABCD",
     [
         (0, 3, 3, 0),  # ABC
         (6, 4, 2, 0),  # ABCD
     ],
     #  0  1  2  3  4  5  6  7  8  9 10 11 12 13
     #  A  B  C  A  B  C  A  B  C  D  A  B  C  D
     #  .  .  .  .  .  0  .  .  0  .  .  .  .  1
     #  0  0  0  0  0  6  6  6  9  9  9  9  9 14
     [[], [], [], [], [], [0], [0], [0], [0], [], [], [], [], [1]],  # flat
     [0, 0, 0, 0, 0, 6, 6, 6, 9, 9, 9, 9, 9, 14],
     {"ABC[2];ABCD[2]"},
     ),

    ("CCACACAC",
     [
         (0, 1, 2, 0),  # C
         (1, 2, 3, 1),  # CA
     ],
     #  0  1  2  3  4  5  6  7
     #  C  C  A  C  A  C  A  C
     #  .  0  .  .  1  1  1  1
     #  0  2  2  2  4  6  6  8
     [[], [0], [], [], [1], [1], [1], [1]],  # flat
     [0, 2, 2, 2, 4, 6, 6, 8],
     {"C[2];AC[3]"},
     ),

    ("CACATCACATCACATCACAT",
     [
         (0, 2, 2, 0),   # CA
         (0, 5, 4, 0),   # CACAT
         (5, 2, 2, 0),   # CA
         (10, 2, 2, 0),  # CA
         (15, 2, 2, 0),  # CA
     ],
     #  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
     #  C  A  C  A  T  C  A  C  A  T  C  A  C  A  T  C  A  C  A  T
     #  .  .  .  0  .  .  .  .  2  1  1  1  1  1  1  1  1  1  1  1
     #                                         3              4
     #  0  0  0  4  4  4  4  4  8 10 10 10 10 14 15 15 15 15 19 20
     [[], [], [], [0], [], [], [], [], [2], [1], [1], [1], [1], [1, 3], [1], [1], [1], [1], [1, 4], [1]],  # nesting
     [0, 0, 0, 4, 4, 4, 4, 4, 8, 10, 10, 10, 10, 14, 15, 15, 15, 15, 19, 20],
     {"CACAT[4]"},
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
     #  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
     #  A  A  C  A  A  C  A  A  A  A  C  A  A  C  A
     #  .  0  .  .  2  1  1  1  4  4  .  .  3  5  5
     #                       4              6
     #  0  2  2  2  4  6  6  8  9 10 10 10 12 14 15
     [[], [0], [], [], [2], [1], [1], [1, 4], [4], [4], [], [], [3, 6], [5], [5]],  # interval intersection & prefix
     [0, 2, 2, 2, 4, 6, 6, 8, 9, 10, 10, 10, 12, 14, 15],
     {"AAC[2];A[3];ACA[2]"},
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
     #  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
     #  A  A  A  A  A  A  C  A  A  A  A  A  A  C  A  A  A  A  A  A  C  A  C  A  C  A  A  A  A  A  A  A
     #  .  0  0  0  0  0  .  .  2  2  2  2  2  1  1  1  1  1  1  1  1  1  4  4  4  4  5  5  5  5  5  5
     #                                               3  3  3  3  3
     #  0  2  3  4  5  6  6  6  8  9 10 11 12 14 14 16 17 18 19 20 21 21 23 24 25 26 27 28 29 30 31 32
     [[], [0], [0], [0], [0], [0], [], [], [2], [2], [2], [2], [2], [1], [1], [1, 3], [1, 3], [1, 3], [1, 3], [1, 3], [1], [1], [4], [4], [4], [4], [5], [5], [5], [5], [5], [5]],  # nesting
     [0, 2, 3, 4, 5, 6, 6, 6, 8, 9, 10, 11, 12, 14, 14, 16, 17, 18, 19, 20, 21, 21, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32],
     {
        "AAAAAAC[2];A[6];CA[3];A[6]",
        "AAAAAAC[2];A[5];AC[3];A[7]",
        "A[5];ACAAAAA[2];AC[3];A[7]",
        "AAAAAAC[3];AC[2];A[7]",
        "A[6];CAAAAAA[2];CA[3];A[6]",
        "A[2];AAAACAA[2];A[4];CA[3];A[6]",
        "A[3];AAACAAA[2];A[3];CA[3];A[6]",
        "A[4];AACAAAA[2];A[2];CA[3];A[6]",
        "A[3];AAACAAA[2];A[2];AC[3];A[7]",
        "A[2];AAAACAA[2];A[3];AC[3];A[7]",
     },
     ),

    ("CATCATACATACTACTAAAAA",
     [
         (0, 3, 2, 0),   # CAT
         (3, 4, 2, 1),   # CATA
         (9, 3, 2, 2),   # TAC
         (16, 1, 5, 0),  # A
     ],
     #  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
     #  C  A  T  C  A  T  A  C  A  T  A  C  T  A  C  T  A  A  A  A  A
     #  .  .  .  .  .  0  .  .  .  .  1  1  .  .  2  2  2  3  3  3  3
     #  0  0  0  0  0  6  6  6  6  6  8  8  8  8 12 12 14 14 16 17 18
     [[], [], [], [], [], [0], [], [], [], [], [1], [1], [], [], [2], [2], [2], [3], [3], [3], [3]],  # flat
     [0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 8, 8, 8, 8, 12, 12, 14, 14, 16, 17, 18],
     {"CAT;CATA[2];CTA[2];A[4]"},
     ),

    ("CTACTACTATATA",
     [
         (0, 3, 3, 0),  # CTA
         (7, 2, 3, 0),  # TA
     ],
     #  0  1  2  3  4  5  6  7  8  9 10 11 12
     #  C  T  A  C  T  A  C  T  A  T  A  T  A
     #  .  .  .  .  .  0  0  0  0  .  1  1  1
     #  0  0  0  0  0  6  6  6  9  9 10 10 13
     [[], [], [], [], [], [0], [0], [0], [0], [], [1], [1], [1]],  # flat
     [0, 0, 0, 0, 0, 6, 6, 6, 9, 9, 10, 10, 13],
     {"CTA[3];TA[2]"},
     ),

    ("AATAAATAAAAA",
     [
         (0, 1, 2, 0),  # A
         (0, 4, 2, 2),  # AATA
         (3, 1, 3, 0),  # A
         (7, 1, 5, 0),  # A
     ],
     #  0  1  2  3  4  5  6  7  8  9 10 11
     #  A  A  T  A  A  A  T  A  A  A  A  A
     #  .  0  .  .  2  2  .  1  1  1  3  3
     #                          3  3
     #  0  2  2  2  4  5  5  8  8 10 11 12
     [[], [0], [], [], [2], [2], [], [1], [1, 3], [1, 3], [3], [3]],  # interval intersection
     [0, 2, 2, 2, 4, 5, 5, 8, 8, 10, 11, 12],
     {
         "AATA[2];A[4]",
         "A[2];TAAA[2];A[2]",
     },
     ),

    ("TAATAAAAAA",
     [
         (0, 3, 2, 0),  # TAA
         (1, 1, 2, 0),  # A
         (4, 1, 6, 0),  # A
     ],
     #  0  1  2  3  4  5  6  7  8  9
     #  T  A  A  T  A  A  A  A  A  A
     #  .  .  1  .  .  0  2  2  2  2
     #                 2
     #  0  0  2  2  2  6  6  8  9 10
     [[], [], [1], [], [], [0, 2], [2], [2], [2], [2]],  # prefix
     [0, 0, 2, 2, 2, 6, 6, 8, 9, 10],
     {"TAA[2];A[4]"},
     ),

    ("CCCCACCAT",
     [
         (0, 1, 4, 0),  # C
         (2, 3, 2, 0),  # CCA
         (5, 1, 2, 0),  # C
     ],
     #  0  1  2  3  4  5  6  7  8
     #  C  C  C  C  A  C  C  A  T
     #  .  0  0  0  .  .  2  1  .
     #  0  2  3  4  4  4  6  8  8
     [[], [0], [0], [0], [], [], [2], [1], []],  # flat
     [0, 2, 3, 4, 4, 4, 6, 8, 8],
     {"C[2];CCA[2];T"},
     ),

    ("TACCACCCC",
     [
         (1, 3, 2, 0),  # ACC
         (2, 1, 2, 0),  # C
         (5, 1, 4, 0),  # C
     ],
     #  0  1  2  3  4  5  6  7  8
     #  T  A  C  C  A  C  C  C  C
     #  .  .  .  1  .  .  0  2  2
     #                    2
     #  0  0  0  2  2  2  6  6  8
     [[], [], [], [1], [], [], [0, 2], [2], [2]],  # prefix
     [0, 0, 0, 2, 2, 2, 6, 6, 8],
     {"T;ACC[2];C[2]"},
     ),

    ("CCACACAAA",
     [
         (0, 1, 2, 0),  # C
         (1, 2, 3, 0),  # CA
         (6, 1, 3, 0),  # A
     ],
     #  0  1  2  3  4  5  6  7  8
     #  C  C  A  C  A  C  A  A  A
     #  .  0  .  .  1  1  1  2  2
     #  0  2  2  2  4  6  6  8  9
     [[], [0], [], [], [1], [1], [1], [2], [2]],  # flat
     [0, 2, 2, 2, 4, 6, 6, 8, 9],
     {"C[2];AC[2];A[3]"},
     ),

    ("CGCGGGGCGC",
     [
         (0, 2, 2, 0),  # CG
         (3, 1, 4, 0),  # G
         (6, 2, 2, 0),  # GC
     ],
     #  0  1  2  3  4  5  6  7  8  9
     #  C  G  C  G  G  G  G  C  G  C
     #  .  .  .  0  1  1  1  .  .  2
     #  0  0  0  4  4  6  7  7  7 10
     [[], [], [], [0], [1], [1], [1], [], [], [2]],  # flat
     [0, 0, 0, 4, 4, 6, 7, 7, 7, 10],
     {"CG[2];G[2];GC[2]"},
     ),

    ("AACAACACAC",
     [
         (0, 1, 2, 0),  # A
         (0, 3, 2, 1),  # AAC
         (3, 1, 2, 0),  # A
         (4, 2, 3, 0),  # AC
     ],
     #  0  1  2  3  4  5  6  7  8  9
     #  A  A  C  A  A  C  A  C  A  C
     #  .  0  .  .  2  1  1  3  3  3
     #  0  2  2  2  4  6  6  6  8 10
     [[], [0], [], [], [2], [1], [1], [3], [3], [3]],  # flat
     [0, 2, 2, 2, 4, 6, 6, 6, 8, 10],
     {"AAC[2];AC[2]"},
     ),

    ("TATAATAATA",
     [
         (0, 2, 2, 0),  # TA
         (1, 3, 3, 0),  # ATA
         (3, 1, 2, 0),  # A
         (6, 1, 2, 0),  # A
     ],
     #  0  1  2  3  4  5  6  7  8  9
     #  T  A  T  A  A  T  A  A  T  A
     #  .  .  .  0  2  .  1  1  1  1
     #                       3
     #  0  0  0  4  4  4  6  6  6 10
     [[], [], [], [0], [2], [], [1], [1, 3], [1], [1]],  # nesting
     [0, 0, 0, 4, 4, 4, 6, 6, 6, 10],
     {"TA[2];ATA[2]"},
     ),

    ("AAACAACACAC",
     [
         (0, 1, 3, 0),  # A
         (4, 1, 2, 0),  # A
         (5, 2, 3, 0),  # AC
         (1, 3, 2, 1),  # AAC
     ],
     #  0  1  2  3  4  5  6  7  8  9 10
     #  A  A  A  C  A  A  C  A  C  A  C
     #  .  0  0  .  .  1  3  3  2  2  2
     #  0  2  3  3  3  5  6  8  8  9 10
     [[], [0], [0], [], [], [1], [3], [3], [2], [2], [2]],  # flat, counter-proof of hwm approach
     [0, 2, 3, 3, 3, 5, 6, 8, 8, 9, 10],
     {"A;AAC[2];AC[2]"},
     ),

    ("CGAAGAAGAC",
     [
        (2, 1, 2, 0),  # A
        (5, 1, 2, 0),  # A
        (1, 3, 2, 2),  # GAA
     ],
     #  0  1  2  3  4  5  6  7  8  9
     #  C  G  A  A  G  A  A  G  A  C
     #  .  .  .  0  .  .  1  2  2  .
     #                    2
     #  0  0  0  2  2  2  6  6  6  6
     [[], [], [], [0], [], [], [1, 2], [2], [2], []],  #
     [0, 0, 0, 2, 2, 2, 6, 6, 6, 6],
     {
         "C;GAA[2];GAC",
         "CGA;AGA[2];C",
         "CG;AAG[2];AC",
     }
    ),

    # ("CACCACACCACCACACCACACCACCACACCACCA",
    #  [
    #     (2, 1, 2, 0),   #  0: C
    #     (7, 1, 2, 0),   #  1: C
    #     (10, 1, 2, 0),  #  2: C
    #     (15, 1, 2, 0),  #  3: C
    #     (20, 1, 2, 0),  #  4: C
    #     (23, 1, 2, 0),  #  5: C
    #     (28, 1, 2, 0),  #  6: C
    #     (31, 1, 2, 0),  #  7: C
    #     (3, 2, 2, 1),   #  8: CA
    #     (11, 2, 2, 1),  #  9: CA
    #     (16, 2, 2, 1),  # 10: CA
    #     (24, 2, 2, 1),  # 11: CA
    #     (0, 3, 2, 0),   # 12: CAC
    #     (5, 3, 3, 0),   # 13: CAC
    #     (13, 3, 2, 0),  # 14: CAC
    #     (18, 3, 3, 0),  # 15: CAC
    #     (26, 3, 2, 2),  # 16: CAC
    #     (0, 5, 2, 1),   # 17: CACCA
    #     (8, 5, 3, 1),   # 18: CACCA
    #     (21, 5, 2, 1),  # 19: CACCA
    #     (0, 8, 2, 3),   # 20: CACCACAC
    #     (13, 8, 2, 5),  # 21: CACCACAC
    #     (0, 13, 2, 6),  # 22: CACCACACCACCA
    #  ],
    #  #  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33
    #  #  C  A  C  C  A  C  A  C  C  A  C  C  A  C  A  C  C  A  C  A  C  C  A  C  C  A  C  A  C  C  A  C  C  A
    #  #  .  .  .  0  . 12  8  8  1 17 13  2 13 13  9  9  3 18 14 10 10  4 18 15  5 15 15 11 11  6 19 16  7 16
    #  #                               17 13          20 20 20 18 18 18 18    18 15 22 22 22 21 21 21 19 16 21
    #  #                                                       20                            22 22 22 21 21
    #  #                                                                                              22
    #  #  0  0  0  2  2  6  6  6  8 10 10 12
    #  [[], [], [], [0], [], [12], [8], [8], [1], [17], [13, 17], [2, 13], [13], [13], [9], [9, 20], [3, 20], [18, 20], [14, 18, 20], [10, 18], [10, 18], [4, 18], [18], [15, 18], [5, 15], [15, 22], [15, 22], [11, 22], [11, 21, 22], [6, 21, 22], [19, 21, 22], [16, 19, 21, 22], [7, 16, 21], [16, 21]],  # complex
    #  [],
    #  {},
    # ),
]


@pytest.mark.parametrize(
"word, pmrs, inv, max_cover, hgvs", TESTS)
def test_extract_repeats(word, pmrs, inv, max_cover, hgvs):
    assert set(pmrs) == set(extract_repeats(word))


@pytest.mark.parametrize(
    "word, pmrs, inv, max_cover, hgvs", TESTS)
def test_repeats_alt(word, pmrs, inv, max_cover, hgvs):
    assert set(pmrs) == set([(start, len(pattern), count, remainder) for pattern, start, count, remainder in repeats_alt(word)])


@pytest.mark.parametrize(
    "word, pmrs, inv, max_cover, hgvs", TESTS)
def test_find_pmrs(word, pmrs, inv, max_cover, hgvs):
    assert set(pmrs) == set(find_pmrs(word))


@pytest.mark.parametrize(
    "word, pmrs, inv, max_cover, hgvs", TESTS)
def test_inv_array(word, pmrs, inv, max_cover, hgvs):
    assert inv == inv_array(len(word), pmrs)


@pytest.mark.parametrize(
    "word, pmrs, inv, max_cover, hgvs", TESTS)
def test_cover(word, pmrs, inv, max_cover, hgvs):
    assert max_cover == cover(word, pmrs)


@pytest.mark.parametrize(
    "word, pmrs, inv, max_cover, hgvs", TESTS)
def test_brute_cover(word, pmrs, inv, max_cover, hgvs):
    assert max_cover[-1] == brute_cover(word, sorted([(start, period, count, remainder) for start, period, count, remainder in pmrs]))


@pytest.mark.parametrize(
    "word, pmrs, inv, max_cover, hgvs", TESTS)
def test_hgvs_solutions(word, pmrs, inv, max_cover, hgvs):
    max_cover, solutions = cover(word, pmrs, hgvs=True)
    assert hgvs == set([path2hgvs(path, word) for path in sol2path(solutions, max_cover)])


@pytest.mark.parametrize(
    "word, pmrs, inv, max_cover, hgvs", TESTS)
def test_hgvs(word, pmrs, inv, max_cover, hgvs):
    overlap = overlapping(pmrs)
    assert hgvs == set([path2hgvs(path, word) for path in sol2path3(inv, pmrs, max_cover, overlap, word)])
