#include <assert.h>     // assert
#include <stddef.h>     // size_t
#include <stdio.h>      // stderr, fprintf
#include <stdlib.h>     // EXIT_SUCCESS
#include <string.h>     // strlen

#include "../include/alloc.h"       // VA_Allocator
#include "../include/std_alloc.h"   // va_std_allocator


#include "../src/array.c"
#include "../src/edit.c"


static void
test_wu_compare(void)
{
    static struct
    {
        char const* const restrict a;
        char const* const restrict b;
        size_t const distance;
    } const tests[] =
    {
        {            "",             "",  0},
        {           "a",            "a",  0},
        {           "a",            "b",  2},
        {        "abcf",        "abcdf",  1},
        {       "abcdf",         "abcf",  1},
        {          "aa",    "aaaaaaaaa",  7},
        {   "aaaaaaaaa",           "aa",  7},
        {  "acbdeacbed", "acebdabbabed",  6},
        {"acebdabbabed",   "acbdeacbed",  6},
        {    "ACCTGACT",    "ATCTTACTT",  5},
        {   "ATCTTACTT",     "ACCTGACT",  5},
    };

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        size_t const m = strlen(tests[i].a);
        size_t const n = strlen(tests[i].b);

        size_t const wu_distance = m > n ?
            wu_compare(va_std_allocator, n, tests[i].b, m, tests[i].a) :
            wu_compare(va_std_allocator, m, tests[i].a, n, tests[i].b);
        assert(wu_distance == tests[i].distance);

        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_wu_compare


int
main(void)
{
    test_wu_compare();
    return EXIT_SUCCESS;
} // main
