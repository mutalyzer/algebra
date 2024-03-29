#include <assert.h>     // assert
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // UINT32_MAX, uint32_t
#include <stdio.h>      // stderr, fprintf
#include <stdlib.h>     // EXIT_SUCCESS


#include "../src/parser.c"


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
    }; // TESTS

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
    }; // TESTS

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
    }; // TESTS

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
    }; // TESTS

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
    }; // TESTS

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
        { 0,       NULL,    "", 0},
        { 0,         "",    "", 0},
        { 1,        "a",    "", 0},
        { 1,        "a",  "aa", 0},
        { 1,        "a",   "a", 1},
        { 3,      "del", "del", 3},
        { 8, "deletion", "del", 3},
    }; // TESTS

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        size_t const tok = match(tests[i].len, tests[i].expression, tests[i].word);
        assert(tok == tests[i].matched);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_match


int
main(void)
{
    test_is_digit();
    test_is_dna_nucleotide();
    test_match_number();
    test_match_sequence();
    test_match_location();
    test_match();
    return EXIT_SUCCESS;
} // main
