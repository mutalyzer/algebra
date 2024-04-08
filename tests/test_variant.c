#include <assert.h>     // assert
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // size_t
#include <stdio.h>      // stderr, fprintf
#include <stdlib.h>     // EXIT_SUCCESS


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


int
main(void)
{
    test_va_variant_eq();
    test_va_variant_len();
    return EXIT_SUCCESS;
} // main
