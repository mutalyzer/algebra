import pytest
from algebra import Variant
from algebra.extractor import to_hgvs
from algebra.extractor.limits import limit_minimal_dfs, limit_minimal_accordion
from algebra.lcs import lcs_graph


@pytest.mark.parametrize(
    "reference, observed, left, right, variants, hgvs",
    [
        (
            "TCCATCCGGAGGGCTCGCC",
            "TCCGATAGCC",
            0,
            17,
            [
                Variant(3, 3, "G"),
                Variant(5, 9, ""),
                Variant(11, 17, ""),
            ],
            "[3_4insG;6_9del;12_17del]",
        ),
        (
                "TCCATCCGGAGGGCTCGCC",
                "TCCGATAGCC",
                1,
                18,
                [
                    Variant(3, 3, "G"),
                    Variant(5, 9, ""),
                    Variant(11, 13, ""),
                    Variant(14, 18, "")
                ],
                "[3_4insG;6_9del;12_13del;15_18del]",
        ),
        (
                "TCCATCCGGAGGGCTCGCC",
                "TCCGATAGCC",
                1,
                17,
                [
                    Variant(3, 3, "G"),
                    Variant(5, 9, ""),
                    Variant(11, 17, ""),
                ],
                "[3_4insG;6_9del;12_17del]",
        ),
        (
                "TCCATCCGGAGGGCTCGCC",
                "TCCGATAGCC",
                2,
                17,
                [
                    Variant(3, 3, "G"),
                    Variant(5, 9, ""),
                    Variant(11, 17, ""),
                ],
                "[3_4insG;6_9del;12_17del]",
        ),
        (
                "TCCATCCGGAGGGCTCGCC",
                "TCCGATAGCC",
                2,
                16,
                [
                    Variant(3, 3, "G"),
                    Variant(5, 9, ""),
                    Variant(10, 16, ""),
                ],
                "[3_4insG;6_9del;11_16del]",
        ),
        (
                "TCCATCCGGAGGGCTCGCC",
                "TCCGATAGCC",
                3,
                16,
                [],
                "=",
        ),
        (
                "TCCATCCGGAGGGCTCGCC",
                "TCCGATAGCC",
                2,
                15,
                [],
                "=",
        ),
        (
            "ATATATTTT",
            "ATTGCCTT",
            0,
            8,
            [
                Variant(2, 3, ""),
                Variant(4, 5, "GCC"),
                Variant(6, 8, "")],
            "[3del;5delinsGCC;7_8del]",
        ),
        (
            "ATATATTTT",
            "ATTGCCTT",
            1,
            8,
            [
                Variant(2, 3, ""),
                Variant(4, 5, "GCC"),
                Variant(6, 8, "")
            ],
            "[3del;5delinsGCC;7_8del]",
        ),
        (
            "ATATATTTT",
            "ATTGCCTT",
            1,
            7,
            [
                Variant(2, 3, ""),
                Variant(4, 7, "GCC"),
            ],
            "[3del;5_7delinsGCC]",
        ),
        (
            "ATATATTTT",
            "ATTGCCTT",
            2,
            7,
            [],
            "=",
        ),
        (
            "ATATATTTT",
            "ATTGCCTT",
            1,
            6,
            [],
            "=",
        ),
    ],
)
def test_restrict_minimal_dfs_first(reference, observed, left, right, variants, hgvs):
    root = lcs_graph(reference, observed)
    minimal = list(limit_minimal_dfs(root, left, right))
    if minimal:
        minimal = minimal[0]
    assert minimal == variants
    assert to_hgvs(minimal, reference) == hgvs



@pytest.mark.parametrize(
    "reference, observed, left, right, variants, hgvs",
    [
        (
            "TCCATCCGGAGGGCTCGCC",
            "TCCGATAGCC",
            0,
            17,
            [
                Variant(3, 8, ""),
                Variant(10, 14, ""),
                Variant(15, 16, "A"),
            ],
            "[4_8del;11_14del;16C>A]",
        ),
        (
                "TCCATCCGGAGGGCTCGCC",
                "TCCGATAGCC",
                1,
                18,
                [
                    Variant(3, 8, ""),
                    Variant(10, 14, ""),
                    Variant(15, 16, "A"),
                ],
                "[4_8del;11_14del;16C>A]",
        ),
        (
                "TCCATCCGGAGGGCTCGCC",
                "TCCGATAGCC",
                1,
                17,
                [
                    Variant(3, 8, ""),
                    Variant(10, 14, ""),
                    Variant(15, 16, "A"),
                ],
                "[4_8del;11_14del;16C>A]",
        ),
        (
                "TCCATCCGGAGGGCTCGCC",
                "TCCGATAGCC",
                2,
                17,
                [
                    Variant(3, 8, ""),
                    Variant(10, 14, ""),
                    Variant(15, 16, "A"),
                ],
                "[4_8del;11_14del;16C>A]",
        ),
        # (
        #         "TCCATCCGGAGGGCTCGCC",
        #         "TCCGATAGCC",
        #         2,
        #         16,
        #         [
        #             Variant(3, 3, "G"),
        #             Variant(5, 9, ""),
        #             Variant(10, 16, ""),
        #         ],
        #         "[3_4insG;6_9del;11_16del]",
        # ),
        (
                "TCCATCCGGAGGGCTCGCC",
                "TCCGATAGCC",
                3,
                16,
                None,
                "=",
        ),
        # (
        #         "TCCATCCGGAGGGCTCGCC",
        #         "TCCGATAGCC",
        #         2,
        #         15,
        #         [],
        #         "=",
        # ),
        # (
        #     "ATATATTTT",
        #     "ATTGCCTT",
        #     0,
        #     8,
        #     [
        #         Variant(2, 3, ""),
        #         Variant(4, 5, "GCC"),
        #         Variant(6, 8, "")],
        #     "[3del;5delinsGCC;7_8del]",
        # ),
        # (
        #     "ATATATTTT",
        #     "ATTGCCTT",
        #     1,
        #     8,
        #     [
        #         Variant(2, 3, ""),
        #         Variant(4, 5, "GCC"),
        #         Variant(6, 8, "")
        #     ],
        #     "[3del;5delinsGCC;7_8del]",
        # ),
        # (
        #     "ATATATTTT",
        #     "ATTGCCTT",
        #     1,
        #     7,
        #     [
        #         Variant(2, 3, ""),
        #         Variant(4, 7, "GCC"),
        #     ],
        #     "[3del;5_7delinsGCC]",
        # ),
    ],
)
def test_restrict_minimal_accordion(reference, observed, left, right, variants, hgvs):
    root = lcs_graph(reference, observed)
    minimal = limit_minimal_accordion(root, left, right)
    assert minimal == variants
    assert to_hgvs(minimal, reference) == hgvs

