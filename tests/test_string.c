#include <assert.h>     // assert
#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // stderr, fprintf
#include <stdlib.h>     // EXIT_SUCCESS
#include <string.h>     // memcmp

#include "../include/alloc.h"       // VA_Allocator
#include "../include/std_alloc.h"   // std_allocator


#include "../src/string.c"


static void
test_va_string_concat(void)
{
    static struct {
        VA_Allocator const* const allocator;
        VA_String lhs;
        VA_String const rhs;
        VA_String const expected;
    } const tests[] =
    {
        {&std_allocator, {0, NULL}, {0, NULL}, {0, NULL}},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        VA_String const result = va_string_concat(tests[i].allocator, tests[i].lhs, tests[i].rhs);
        assert(memcmp(result.data, tests[i].expected.data, tests[i].expected.length) == 0);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_string_concat


static void
test_va_string_slice(void)
{
    static struct {
        VA_Allocator const* const allocator;
        VA_String const string;
        size_t const start;
        size_t const end;
        VA_String const expected;
    } const tests[] =
    {
        {&std_allocator, {0, NULL}, 0, 0, {0, NULL}},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        VA_String const result = va_string_slice(tests[i].allocator, tests[i].string, tests[i].start, tests[i].end);
        assert(memcmp(result.data, tests[i].expected.data, tests[i].expected.length) == 0);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_string_slice


int
main(void)
{
    test_va_string_concat();
    test_va_string_slice();
    return EXIT_SUCCESS;
} // main
