#include <assert.h>     // assert
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // UINT32_MAX, uint32_t
#include <stdio.h>      // stderr, fprintf
#include <stdlib.h>     // EXIT_SUCCESS
#include <string.h>     // memcmp


#include "../src/parser.c"


#define MAX_VARIANTS 2


static void
test_is_digit(void)
{
    static struct {
        char const ch;
        bool const expected;
    } const tests[] =
    {
        {  '\0', false},
        {   '/', false},
        {   '0',  true},
        {   '1',  true},
        {   '2',  true},
        {   '3',  true},
        {   '4',  true},
        {   '5',  true},
        {   '6',  true},
        {   '7',  true},
        {   '8',  true},
        {   '9',  true},
        {   ':', false},
        {   'a', false},
        {   'A', false},
        {'\xFF', false},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        assert(tests[i].expected == is_digit(tests[i].ch));
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_is_digit


static void
test_is_dna_nucleotide(void)
{
    static struct {
        char const ch;
        bool const expected;
    } const tests[] =
    {
        {  '\0', false},
        {   '0', false},
        {   'a', false},
        {   'A',  true},
        {   'C',  true},
        {   'G',  true},
        {   'T',  true},
        {   'B', false},
        {'\xFF', false},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        assert(tests[i].expected == is_dna_nucleotide(tests[i].ch));
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_is_dna_nucleotide


static void
test_match_number(void)
{
    static struct {
        size_t const      len;
        char const* const expression;
        size_t const      matched;
        uint32_t const    expected;
    } const tests[] =
    {
        { 0,         NULL,  0,          0},
        { 0,           "",  0,          0},
        { 1,          "1",  1,          1},
        {10,         "10",  2,         10},  // saved by \0
        { 1,         "10",  1,          1},
        {10, "123a123a12",  3,        123},
        { 5,      "-1000",  0,          0},
        {10, "1234567890", 10, 1234567890},
        {10, "4294967295", 10, UINT32_MAX},
        {10, "4294967296", 10,          0},  // overflow
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        uint32_t num = 0;
        size_t const tok = match_number(tests[i].len, tests[i].expression, &num);
        assert(tok == tests[i].matched);
        assert(num == tests[i].expected);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_match_number


static void
test_match_sequence(void)
{
    static struct {
        size_t const      len;
        char const* const expression;
        size_t const      matched;
    } const tests[] =
    {
        {0,   NULL, 0},
        {0,     "", 0},
        {1,    "A", 1},
        {4, "ACGT", 4},
        {4, "acgt", 0},
        {4, "ACgt", 2},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        size_t const tok = match_sequence(tests[i].len, tests[i].expression);
        assert(tok == tests[i].matched);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_match_sequence


static void
test_match_location(void)
{
    static struct {
        size_t const      len;
        char const* const expression;
        size_t const      matched;
        uint32_t const    expected_start;
        uint32_t const    expected_end;
    } const tests[] =
    {
        { 0,         NULL,  0,    0,     0},
        { 0,           "",  0,    0,     0},
        { 1,          "0",  1,    0,     0},
        { 1,          "1",  1,    1,     1},
        { 1,          "a",  0,    0,     0},
        { 2,         "9a",  1,    9,     9},
        { 3,        "0_1",  3,    0,     1},
        { 2,         "0_",  0,    0,     0},
        { 3,        "0_a",  0,    0,     0},
        { 3,        "1_1",  0,    1,     1},
        { 3,        "1_0",  0,    1,     0},
        {10, "1234_56789", 10, 1234, 56789},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        uint32_t start = 0;
        uint32_t end = 0;
        size_t const tok = match_location(tests[i].len, tests[i].expression, &start, &end);
        assert(tok == tests[i].matched);
        assert(start == tests[i].expected_start);
        assert(end == tests[i].expected_end);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_match_location


static void
test_match(void)
{
    static struct {
        size_t const      len;
        char const* const expression;
        char const* const word;
        size_t const      matched;
    } const tests[] =
    {
        {0,       NULL,         "", 0},
        {0,         "",         "", 0},
        {1,        "a",         "", 0},
        {1,        "a",       "aa", 0},
        {1,        "a",        "a", 1},
        {3,      "del",      "del", 3},
        {8, "deletion",      "del", 3},
        {3,      "del", "deletion", 0},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        size_t const tok = match(tests[i].len, tests[i].expression, tests[i].word);
        assert(tok == tests[i].matched);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_match


static void
test_match_until(void)
{
    static struct {
        size_t const      len;
        char const* const expression;
        char const        delim;
        size_t const      matched;
    } const tests[] =
    {
        {0,       NULL, '\0', 0},
        {0,         "", '\0', 0},
        {1,        "1",  ':', 0},
        {1,        ":",  ':', 0},
        {8, "AATT:ATT",  ':', 4},
        {3,      ":::",  ':', 0},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        size_t const tok = match_until(tests[i].len, tests[i].expression, tests[i].delim);
        assert(tok == tests[i].matched);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_match_until


static void
test_match_variant(void)
{
    static struct {
        size_t const      len;
        char const* const expression;
        size_t const      matched;
        VA_Variant const  variant;
    } const tests[] =
    {
        {0,        NULL, 0, {0}},
        {0,          "", 0, {0}},
        {1,         "1", 0, {0, 1, 0, 0}},
        {4,      "1del", 4, {0, 1, 0, 0}},
        {4,      "0del", 0, {0}},
        {5,     "1delT", 5, {0, 1, 0, 0}},
        {7,   "1delins", 0, {0, 1, 7, 0}},
        {8,  "1delinsA", 8, {0, 1, 7, 8}},
        {9, "1delAinsA", 9, {0, 1, 8, 9}},
        {6,    "1_2ins", 0, {1, 1, 6, 0}},
        {4,      "1ins", 0, {1, 1, 0, 0}},
        {6,    "0_1ins", 0, {0, 0, 6, 0}},
        {7,   "0_1insA", 7, {0, 0, 6, 7}},
        {6,    "1_3ins", 0, {1, 3, 0, 0}},
        {7,   "1_2insT", 7, {1, 1, 6, 7}},
        {2,        "0>", 0, {0}},
        {2,        "1>", 0, {0, 1, 2, 0}},
        {3,       "1A>", 0, {0, 1, 3, 0}},
        {4,      "1A>T", 4, {0, 1, 3, 4}},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        VA_Variant variant = {0};
        size_t const tok = match_variant(tests[i].len, tests[i].expression, &variant);
        assert(tok == tests[i].matched);
        assert(memcmp(&variant, &tests[i].variant, sizeof(variant)) == 0);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_match_variant


static void
test_va_parse_hgvs(void)
{
    static struct {
        size_t const      len;
        char const* const expression;
        size_t const      count;
        VA_Variant const  variants[MAX_VARIANTS];
    } const tests[] =
    {
        { 0,               NULL, 0, {{0}}},
        { 0,                 "", 0, {{0}}},
        { 1,                "A", 0, {{0}}},
        { 1,                ":", 0, {{0}}},
        { 1,                "g", 0, {{0}}},
        { 2,               "g.", 0, {{0}}},
        { 3,              ":g.", 0, {{0}}},
        { 1,                "1", 0, {{0}}},
        { 3,              "A:1", 0, {{0}}},
        { 5,            "A:g.1", 0, {{0}}},
        { 7,          "A:g.1_2", 0, {{0}}},
        { 1,                "=", 1, {{0}}},
        { 3,              "g.=", 1, {{0}}},
        { 4,             "1del", 1, {{0, 1, 0, 0}}},
        { 6,           "g.1del", 1, {{0, 1, 0, 0}}},
        { 6,           "A:1del", 1, {{0, 1, 0, 0}}},
        { 8,         "A:g.1del", 1, {{0, 1, 0, 0}}},
        { 5,            "1del ", 0, {{0}}},
        { 1,                "[", 0, {{0}}},
        { 5,            "[1del", 0, {{0}}},
        { 6,           "[1del;", 0, {{0}}},
        { 6,           "[1del]", 1, {{0, 1, 0, 0}}},
        { 7,          "[1del] ", 0, {{0}}},
        {10,       "[1del;2del", 0, {{0}}},
        {11,      "[1del;2del]", 2, {{0, 1, 0, 0}, {1, 2, 0, 0}}},
        {10,       "[1del2del]", 0, {{0}}},
        {16, "[1del;2del;3del]", 0, {{0}}},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        VA_Variant variants[MAX_VARIANTS] = {{0}};
        size_t const tok = va_parse_hgvs(tests[i].len, tests[i].expression, MAX_VARIANTS, variants);
        assert(tok == tests[i].count);
        assert(memcmp(&variants, &tests[i].variants, sizeof(variants[0]) * tests[i].count) == 0);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_parse_hgvs


static void
test_va_parse_spdi(void)
{
    static struct {
        size_t const      len;
        char const* const expression;
        size_t const      count;
        VA_Variant const  variants;
    } const tests[] =
    {
        {0,       NULL, 0, {0}},
        {0,         "", 0, {0}},
        {1,        "A", 0, {0}},
        {1,        ":", 0, {0}},
        {2,       "A:", 0, {0}},
        {3,      "A:1", 0, {1, 0, 0, 0}},
        {4,     "A:1:", 0, {1, 0, 0, 0}},
        {5,    "A:1:1", 0, {1, 2, 0, 0}},
        {5,    "A:1:A", 0, {1, 2, 0, 0}},
        {6,   "A:1:1:", 0, {1, 2, 6, 0}},
        {6,   "A:1:A:", 0, {1, 2, 6, 0}},
        {7,  "A:1:1:A", 1, {1, 2, 6, 7}},
        {7,  "A:1:A:A", 1, {1, 2, 6, 7}},
        {8, "A:1:A:A ", 0, {1, 2, 6, 7}},
        {6,   ":1:A:A", 1, {1, 2, 5, 6}},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        VA_Variant variants[1] = {{0}};
        size_t const tok = va_parse_spdi(tests[i].len, tests[i].expression, variants);
        assert(tok == tests[i].count);
        assert(memcmp(&variants, &tests[i].variants, sizeof(variants[0])) == 0);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_parse_spdi


int
main(void)
{
    test_is_digit();
    test_is_dna_nucleotide();
    test_match_number();
    test_match_sequence();
    test_match_location();
    test_match();
    test_match_until();
    test_match_variant();
    test_va_parse_hgvs();
    test_va_parse_spdi();
    return EXIT_SUCCESS;
} // main
