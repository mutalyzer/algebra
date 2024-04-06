#include <assert.h>     // assert
#include <stddef.h>     // size_t
#include <stdio.h>      // stderr, fprintf
#include <stdlib.h>     // EXIT_SUCCESS


#include "../src/variant.c"


static void
test_va_variant_len(void)
{
    static struct {
        VA_Variant const variant;
        size_t const     len;
    } const tests[] =
    {
        {{0},                                          0},
        {{.start=0, .end=0, .obs_start=0, .obs_end=0}, 0},
        {{.start=0, .end=1, .obs_start=0, .obs_end=0}, 1},
        {{.start=0, .end=0, .obs_start=0, .obs_end=1}, 1},
        {{.start=0, .end=1, .obs_start=0, .obs_end=1}, 2},
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
    test_va_variant_len();
    return EXIT_SUCCESS;
} // main
