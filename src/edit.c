#include <inttypes.h>   // intmax_t
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t


#include "../include/allocator.h"   // GVA_Allocator
#include "../include/edit.h"        // gva_edit_distance
#include "common.h"     // ABS, MAX


//#include <stdio.h>  // DEBUG


size_t max_lcs_pos = 0;
MNode* matches = NULL;
char* uniq = NULL;


static inline size_t
onp_snake(size_t const m, char const a[static restrict m],
    size_t const n, char const b[static restrict n],
    intmax_t const k,
    intmax_t const lower,
    intmax_t const upper,
    intmax_t const delta,
    size_t const p,
    bool const swapped)
{
    size_t col = MAX(lower, upper);
    size_t row = col - k;

    intmax_t const d_row = m - row;
    intmax_t const d_col = n - col;

    while (row < m && col < n && a[row] == b[col])
    {
        row += 1;
        col += 1;
        size_t const lcs_pos = (row + col - ABS(delta) - 2 * p + ABS(d_row - d_col)) / 2;
        //fprintf(stderr, "%zu: (%zu, %zu)\n", lcs_pos - 1, row - 1, col - 1);
        uniq[lcs_pos - 1] = 0;
        //matches[lcs_pos - 1] = (MNode) {-1, -1};
        if (lcs_pos > max_lcs_pos)
        {
            max_lcs_pos = lcs_pos;
            uniq[lcs_pos - 1] = 1;
            if (swapped)
            {
                matches[lcs_pos - 1] = (MNode) {col - 1, row - 1};
            } // if
            else
            {
                matches[lcs_pos - 1] = (MNode) {row - 1, col - 1};
            } // else
        } // if
    } // while
    return col;
} // onp_snake


static size_t
onp_compare(GVA_Allocator const allocator,
    size_t const m, char const a[static restrict m],
    size_t const n, char const b[static restrict n],
    bool const swapped)
{
    intmax_t const delta = n - m;
    size_t const offset = m + 1;
    size_t const size = m + n + 3;
    intmax_t* const restrict fp = allocator.allocate(allocator.context, NULL, 0, size * sizeof(*fp));
    if (fp == NULL)
    {
        return -1;
    } // if

    for (intmax_t i = -1; i <= delta; ++i)
    {
        fp[i + offset] = -1;
    } // for

    max_lcs_pos = 0;
    matches = allocator.allocate(allocator.context, NULL, 0, sizeof(*matches) * m);
    uniq = allocator.allocate(allocator.context, NULL, 0, sizeof(*uniq) * m);

    size_t p = 0;
    while ((size_t) fp[delta + offset] != n)
    {
        fp[-p - 1 + offset] = -1;
        fp[delta + p + 1 + offset] = -1;
        for (intmax_t k = -p; k <= delta - 1; ++k)
        {
            fp[k + offset] = onp_snake(m, a, n, b, k, fp[k - 1 + offset] + 1, fp[k + 1 + offset], delta, p, swapped);
        } // for
        for (intmax_t k = delta + p; k >= delta + 1; --k)
        {
            fp[k + offset] = onp_snake(m, a, n, b, k, fp[k - 1 + offset] + 1, fp[k + 1 + offset], delta, p, swapped);
        } // for
        fp[delta + offset] = onp_snake(m, a, n, b, delta, fp[delta - 1 + offset] + 1, fp[delta + 1 + offset], delta, p, swapped);
        p += 1;
    } // while

    allocator.allocate(allocator.context, fp, size * sizeof(*fp), 0);
    return delta + 2 * (p - 1);
} // onp_compare


inline size_t
gva_edit_distance(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs])
{
    //fprintf(stderr, "ALIGN\n");
    return len_ref > len_obs ?
           onp_compare(allocator, len_obs, observed, len_ref, reference, true) :
           onp_compare(allocator, len_ref, reference, len_obs, observed, false);
} // gva_edit_distance
