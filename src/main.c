#include <assert.h>     // assert
#include <inttypes.h>   // imaxabs
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, ptrdiff_t, size_t
#include <stdio.h>      // stderr, fprintf, printf
#include <stdlib.h>     // EXIT_*
#include <string.h>     // strlen

#include "../include/alloc.h"           // VA_Allocator
#include "../include/array.h"           // va_array_*
#include "../include/static_alloc.h"    // va_static_allocator
#include "../include/std_alloc.h"       // va_std_allocator


static size_t count = 0;


static inline size_t
wu_snake(size_t const m,
         char const a[static restrict m],
         size_t const n,
         char const b[static restrict n],
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
wu_compare(VA_Allocator const allocator[static restrict 1],
           size_t const m,
           char const a[static restrict m],
           size_t const n,
           char const b[static restrict n])
{
    ptrdiff_t const delta = n - m;
    size_t const offset = m + 1;
    size_t const size = m + n + 3;
    ptrdiff_t* const restrict fp = allocator->alloc(allocator->context, NULL, 0, size * sizeof(*fp));
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


typedef struct
{
    size_t row;
    size_t col;
    size_t length;
} LCS_Node;


static inline size_t
umax(size_t const a, size_t const b)
{
    return a > b ? a : b;
} // umax


static inline size_t
umin(size_t const a, size_t const b)
{
    return a < b ? a : b;
} // umin


typedef struct
{
    VA_Allocator const* const restrict allocator;
    size_t const len_ref;
    char const* const restrict reference;
    size_t const len_obs;
    char const* const restrict observed;
    size_t* const restrict diagonals;
    LCS_Node** const restrict lcs_nodes;
} Context;


static size_t
expand(Context const context,
       ptrdiff_t const idx,
       size_t const p,
       size_t* const max_lcs_pos)
{
    ptrdiff_t const delta = context.len_obs - context.len_ref;
    size_t const offset = context.len_ref + 1;

    size_t row;
    size_t col;
    size_t end;
    if (idx > 0)
    {
        row = context.diagonals[offset + idx];
        col = row + idx;
        end = umax(context.diagonals[offset + idx - 1] - 1, context.diagonals[offset + idx + 1]);
    } // if
    else if (idx < 0)
    {
        col = context.diagonals[offset + idx];
        row = col - idx;
        end = umax(context.diagonals[offset + idx - 1], context.diagonals[offset + idx + 1] - 1);
    } // if
    else
    {
        row = context.diagonals[offset + idx];
        col = row + idx;
        end = umax(context.diagonals[offset + idx - 1], context.diagonals[offset + idx + 1]);
    } // else

    size_t steps = end + 1;

    bool matching = false;
    size_t match_row = 0;
    size_t match_col = 0;

    for (size_t i = context.diagonals[offset + idx]; i < end; ++i)
    {
        if (context.reference[row] == context.observed[col])
        {
            if (!matching)
            {
                match_row = row;
                match_col = col;
                matching = true;
            } // if
        } // if
        else if (matching)
        {
            ptrdiff_t const d_row = context.len_ref - row;
            ptrdiff_t const d_col = context.len_obs - col;
            size_t const lcs_pos = (row + col - imaxabs(delta) - 2 * p + imaxabs(d_row - d_col)) / 2 - 1;
            if (lcs_pos > *max_lcs_pos)
            {
                *max_lcs_pos = lcs_pos + 1;
            } // if
            size_t const length = row - match_row;
            if (context.lcs_nodes[lcs_pos] == NULL)
            {
                context.lcs_nodes[lcs_pos] = va_array_init(context.allocator, 1, sizeof(*context.lcs_nodes[lcs_pos]));
            } // if
            va_array_append(context.allocator, context.lcs_nodes[lcs_pos], ((LCS_Node) {match_row, match_col, length}));
            matching = false;
        } // if
        row += 1;
        col += 1;
    } // for

    if (!matching)
    {
        match_row = row;
        match_col = col;
    } // if
    while (row < context.len_ref && col < context.len_obs && context.reference[row] == context.observed[col])
    {
        matching = true;
        row += 1;
        col += 1;
        steps += 1;
    } // while
    if (matching)
    {
        ptrdiff_t const d_row = context.len_ref - row;
        ptrdiff_t const d_col = context.len_obs - col;
        size_t const lcs_pos = (row + col - imaxabs(delta) - 2 * p + imaxabs(d_row - d_col)) / 2 - 1;
        if (lcs_pos > *max_lcs_pos)
        {
            *max_lcs_pos = lcs_pos + 1;
        } // if
        size_t const length = row - match_row;
        if (context.lcs_nodes[lcs_pos] == NULL)
        {
            context.lcs_nodes[lcs_pos] = va_array_init(context.allocator, 1, sizeof(*context.lcs_nodes[lcs_pos]));
        } // if
        va_array_append(context.allocator, context.lcs_nodes[lcs_pos], ((LCS_Node) {match_row, match_col, length}));
    } // if

    return steps;
} // expand


static size_t
edit(VA_Allocator const allocator[static restrict 1],
     size_t const len_ref,
     char const reference[static restrict len_ref],
     size_t const len_obs,
     char const observed[static restrict len_obs])
{
    ptrdiff_t const delta = len_obs - len_ref;
    size_t const offset = len_ref + 1;
    size_t const size = len_ref + len_obs + 3;

    Context const context =
    {
        allocator,
        len_ref,
        reference,
        len_obs,
        observed,
        .diagonals = allocator->alloc(allocator->context, NULL, 0, size * sizeof(*context.diagonals)),
        .lcs_nodes = allocator->alloc(allocator->context, NULL, 0, umin(len_ref, len_obs) * sizeof(*context.lcs_nodes)),
    };
    if (context.diagonals == NULL || context.lcs_nodes == NULL)
    {
        allocator->alloc(allocator->context, context.diagonals, size * sizeof(*context.diagonals), 0);
        allocator->alloc(allocator->context, context.lcs_nodes, umin(len_ref, len_obs) * sizeof(*context.lcs_nodes), 0);
        return -1;
    } // if

    size_t const lower = delta > 0 ? 0 : delta;
    size_t const upper = delta > 0 ? delta : 0;

    size_t const len = umax(len_ref, len_obs) - imaxabs(delta);
    size_t p = 0;
    size_t max_lcs_pos = 0;
    while (context.diagonals[offset + delta] <= len)
    {
        for (ptrdiff_t idx = lower - p; idx < delta; ++idx)
        {
            context.diagonals[offset + idx] = expand(context, idx, p, &max_lcs_pos);
        } // for
        for (ptrdiff_t idx = upper + p; idx > delta; --idx)
        {
            context.diagonals[offset + idx] = expand(context, idx, p, &max_lcs_pos);
        } // for
        context.diagonals[offset + delta] = expand(context, delta, p, &max_lcs_pos);

        p += 1;
    } // while

    allocator->alloc(allocator->context, context.diagonals, size * sizeof(*context.diagonals), 0);

    for (size_t i = 0; i < max_lcs_pos; ++i)
    {
        if (context.lcs_nodes[i] != NULL)
        {
            printf("%zu: ", i);
            for (size_t j = 0; j < va_array_length(context.lcs_nodes[i]); ++j)
            {
                printf("(%zu, %zu, %zu) ", context.lcs_nodes[i][j].row, context.lcs_nodes[i][j].col, context.lcs_nodes[i][j].length);
            } // for
            printf("\n");
        } // if
        context.lcs_nodes[i] = va_array_destroy(allocator, context.lcs_nodes[i]);
    } // for

    allocator->alloc(allocator->context, context.lcs_nodes, umin(len_ref, len_obs) * sizeof(*context.lcs_nodes), 0);

    return imaxabs(delta) + 2 * p - 2;
} // edit


static void
test_wu_compare(void)
{
    static struct
    {
        char const* const restrict a;
        char const* const restrict b;
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

        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_wu_compare


int
main(int argc, char* argv[argc + 1])
{
    // test_wu_compare();

    if (argc < 3)
    {
        fprintf(stderr, "usage: %s reference observed\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    size_t const m = strlen(argv[1]);
    size_t const n = strlen(argv[2]);

    count = 0;
    size_t const edit_distance = edit(&va_std_allocator, m, argv[1], n, argv[2]);
    printf("edit distance: %zu (%zu)\n", edit_distance, count);

    return EXIT_SUCCESS;
} // main
