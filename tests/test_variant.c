#include <assert.h>     // assert
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // size_t
#include <stdio.h>      // stderr, fprintf
#include <stdlib.h>     // EXIT_SUCCESS
#include <string.h>     // strncmp


#include "../src/variant.c"


static void
test_va_variant_eq(void)
{
    static struct {
        VA_Variant const lhs;
        VA_Variant const rhs;
        bool const       eq;
    } const tests[] =
    {
        {         {0},          {0},  true},
        {{0, 0, 0, 0}, {0, 0, 0, 0},  true},
        {{0, 0, 0, 0}, {0, 0, 0, 1}, false},
        {{0, 0, 0, 0}, {0, 0, 1, 0}, false},
        {{0, 0, 0, 0}, {0, 1, 0, 0}, false},
        {{0, 0, 0, 0}, {1, 0, 0, 0}, false},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        bool const eq = va_variant_eq(tests[i].lhs, tests[i].rhs);
        assert(eq == tests[i].eq);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_variant_eq


static void
test_va_variant_len(void)
{
    static struct {
        VA_Variant const variant;
        size_t const     len;
    } const tests[] =
    {
        {         {0}, 0},
        {{0, 0, 0, 0}, 0},
        {{0, 1, 0, 0}, 1},
        {{0, 0, 0, 1}, 1},
        {{0, 1, 0, 1}, 2},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        size_t const len = va_variant_len(tests[i].variant);
        assert(len == tests[i].len);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_variant_len


static void
test_va_patch(void)
{
    static struct {
        size_t const      len_ref;
        char const* const reference;
        size_t const      len_obs;
        char const* const observed;
        size_t const      len_var;
        VA_Variant const  variants[256];
        char const* const expected;
    } const tests[] =
    {
        {5, "AAAAA", 0, "", 0, {{0}}, "AAAAA"},
        {5, "AAAAA", 0, "", 1, {{0, 1, 0, 0}}, "AAAA"},
        {5, "AAAAA", 1, "T", 1, {{0, 1, 0, 1}}, "TAAAA"},
        {5, "AAAAA", 1, "T", 1, {{2, 3, 0, 1}}, "AATAA"},
        {5, "AAAAA", 1, "T", 1, {{8, 9, 0, 1}}, ""},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        char buffer[5] = {'\0'};
        size_t const len = va_patch(tests[i].len_ref, tests[i].reference,
                                    tests[i].len_obs, tests[i].observed,
                                    tests[i].len_var, tests[i].variants,
                                    5, buffer);
        //printf("%zu\n", len);
        assert(strncmp(buffer, tests[i].expected, len) == 0);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_patch


int
main(void)
{
    test_va_variant_eq();
    test_va_variant_len();
    test_va_patch();
    return EXIT_SUCCESS;
} // main
