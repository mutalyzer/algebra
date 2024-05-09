#include <assert.h>     // assert
#include <inttypes.h>   // imaxabs
#include <stddef.h>     // NULL, ptrdiff_t, size_t
#include <stdio.h>      // stderr, fprintf, printf
#include <stdlib.h>     // EXIT_*
#include <string.h>     // strlen

#include "../include/alloc.h"           // VA_Allocator
#include "../include/static_alloc.h"    // va_static_allocator


static size_t count = 0;


static inline size_t
wu_snake(size_t const m,
      char const a[static m],
      size_t const n,
      char const b[static n],
      ptrdiff_t const k,
      ptrdiff_t const lower,
      ptrdiff_t const upper)
{
    size_t col = lower > upper ? lower : upper;
    size_t row = col - k;
    while (row < m && col < n && a[row] == b[col])
    {
        count += 1;
        row += 1;
        col += 1;
    } // while
    if (row < m && col < n)
    {
        count += 1;
    } // if
    return col;
} // wu_snake


static size_t
wu_compare(VA_Allocator const allocator[static 1],
           size_t const m,
           char const a[static m],
           size_t const n,
           char const b[static n])
{
    ptrdiff_t const delta = n - m;
    size_t const offset = m + 1;
    size_t const size = m + n + 3;
    ptrdiff_t* const fp = allocator->alloc(allocator->context, NULL, 0, size * sizeof(*fp));
    if (fp == NULL)
    {
        return -1;
    } // if

    for (size_t i = 0; i < size; ++i)
    {
        fp[i] = -1;
    } // for

    size_t p = 0;
    while ((size_t) fp[delta + offset] != n)
    {
        for (ptrdiff_t k = -p; k <= delta - 1; ++k)
        {
            fp[k + offset] = wu_snake(m, a, n, b, k, fp[k - 1 + offset] + 1, fp[k + 1 + offset]);
        } // for
        for (ptrdiff_t k = delta + p; k >= delta + 1; --k)
        {
            fp[k + offset] = wu_snake(m, a, n, b, k, fp[k - 1 + offset] + 1, fp[k + 1 + offset]);
        } // for
        fp[delta + offset] = wu_snake(m, a, n, b, delta, fp[delta - 1 + offset] + 1, fp[delta + 1 + offset]);
        p += 1;
    } // while

    allocator->alloc(allocator->context, fp, size * sizeof(*fp), 0);

    return delta + 2 * (p - 1);
} // wu_compare


static inline size_t
umax(size_t const a, size_t const b)
{
    return a > b ? a : b;
} // umax


static size_t
expand(size_t const len_ref,
       char const reference[static len_ref],
       size_t const len_obs,
       char const observed[static len_obs],
       ptrdiff_t const idx,
       size_t const offset,
       size_t const diagonals[static len_ref + len_obs + 3])
{
    size_t row;
    size_t col;
    size_t end;
    if (idx > 0)
    {
        row = diagonals[offset + idx];
        col = row + idx;
        end = umax(diagonals[offset + idx - 1] - 1, diagonals[offset + idx + 1]);
    } // if
    else if (idx < 0)
    {
        col = diagonals[offset + idx];
        row = col - idx;
        end = umax(diagonals[offset + idx - 1], diagonals[offset + idx + 1] - 1);
    } // if
    else
    {
        row = diagonals[offset + idx];
        col = row + idx;
        end = umax(diagonals[offset + idx - 1], diagonals[offset + idx + 1]);
    } // else

    row += end - diagonals[offset + idx];
    col += end - diagonals[offset + idx];

    size_t steps = end + 1;
    while (row < len_ref && col < len_obs && reference[row] == observed[col])
    {
        count += 1;
        row += 1;
        col += 1;
        steps += 1;
    } // while
    if (row < len_ref && col < len_obs)
    {
        count += 1;
    } // if

    return steps;
} // expand


static size_t
edit(VA_Allocator const allocator[static 1],
     size_t const len_ref,
     char const reference[static len_ref],
     size_t const len_obs,
     char const observed[static len_obs])
{
    ptrdiff_t const delta = len_obs - len_ref;
    size_t const offset = len_ref + 1;
    size_t const size = len_ref + len_obs + 3;

    size_t* const diagonals = allocator->alloc(allocator->context, NULL, 0, size * sizeof(*diagonals));
    if (diagonals == NULL)
    {
        return -1;
    } // if

    size_t const lower = delta > 0 ? 0 : delta;
    size_t const upper = delta > 0 ? delta : 0;

    size_t const len = umax(len_ref, len_obs) - imaxabs(delta);
    size_t it = 0;
    while (diagonals[offset + delta] <= len)
    {
        for (ptrdiff_t idx = lower - it; idx < delta; ++idx)
        {
            diagonals[offset + idx] = expand(len_ref, reference, len_obs, observed, idx, offset, diagonals);
        } // for
        for (ptrdiff_t idx = upper + it; idx > delta; --idx)
        {
            diagonals[offset + idx] = expand(len_ref, reference, len_obs, observed, idx, offset, diagonals);
        } // for
        diagonals[offset + delta] = expand(len_ref, reference, len_obs, observed, delta, offset, diagonals);

        it += 1;
    } // while

    allocator->alloc(allocator->context, diagonals, size * sizeof(*diagonals), 0);

    return imaxabs(delta) + 2 * it - 2;
} // edit


static void
test_wu_compare(void)
{
    static struct
    {
        char const* const a;
        char const* const b;
        size_t const distance;
        size_t const count;
    } const tests[] =
    {
        {            "",             "",  0,  0},
        {           "a",            "a",  0,  1},
        {           "a",            "b",  2,  1},
        {        "abcf",        "abcdf",  1,  5},
        {       "abcdf",         "abcf",  1,  5},
        {          "aa",    "aaaaaaaaa",  7,  2},
        {   "aaaaaaaaa",           "aa",  7,  2},
        {  "acbdeacbed", "acebdabbabed",  6, 23},
        {"acebdabbabed",   "acbdeacbed",  6, 23},
        {    "ACCTGACT",    "ATCTTACTT",  5, 17},
        {   "ATCTTACTT",     "ACCTGACT",  5, 17},
    };

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        size_t const m = strlen(tests[i].a);
        size_t const n = strlen(tests[i].b);

        count = 0;
        size_t const wu_distance = m > n ?
            wu_compare(&va_static_allocator, n, tests[i].b, m, tests[i].a) :
            wu_compare(&va_static_allocator, m, tests[i].a, n, tests[i].b);
        assert(wu_distance == tests[i].distance);
        assert(count == tests[i].count);

        count = 0;
        size_t const edit_distance = edit(&va_static_allocator, m, tests[i].a, n, tests[i].b);
        assert(edit_distance == tests[i].distance);
        assert(count == tests[i].count);

        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_wu_compare


int
main(int argc, char* argv[argc + 1])
{
    test_wu_compare();

    if (argc < 3)
    {
        fprintf(stderr, "usage: %s reference observed\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    size_t const m = strlen(argv[1]);
    size_t const n = strlen(argv[2]);

    count = 0;
    size_t const wu_distance = m > n ?
        wu_compare(&va_static_allocator, n, argv[2], m, argv[1]) :
        wu_compare(&va_static_allocator, m, argv[1], n, argv[2]);
    printf("  wu distance: %zu (%zu)\n", wu_distance, count);

    count = 0;
    size_t const edit_distance = edit(&va_static_allocator, m, argv[1], n, argv[2]);
    printf("edit distance: %zu (%zu)\n", edit_distance, count);

    return EXIT_SUCCESS;
} // main
