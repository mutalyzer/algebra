#include <inttypes.h>   // intmax_t
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t


#include "../include/allocator.h"   // GVA_Allocator
#include "../include/edit.h"        // gva_edit_distance
#include "common.h"     // ABS, MAX


static inline size_t
onp_snake(size_t const m, char const a[static restrict m],
    size_t const n, char const b[static restrict n],
    intmax_t const k, intmax_t const lower, intmax_t const upper)
{
    size_t col = MAX(lower, upper);
    size_t row = col - k;

    while (row < m && col < n && a[row] == b[col])
    {
        row += 1;
        col += 1;
    } // while
    return col;
} // onp_snake


static size_t
onp_compare(GVA_Allocator const allocator,
    size_t const m, char const a[static restrict m],
    size_t const n, char const b[static restrict n])
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

    size_t p = 0;
    while ((size_t) fp[delta + offset] != n)
    {
        fp[-p - 1 + offset] = -1;
        fp[delta + p + 1 + offset] = -1;
        for (intmax_t k = -p; k <= delta - 1; ++k)
        {
            fp[k + offset] = onp_snake(m, a, n, b, k, fp[k - 1 + offset] + 1, fp[k + 1 + offset]);
        } // for
        for (intmax_t k = delta + p; k >= delta + 1; --k)
        {
            fp[k + offset] = onp_snake(m, a, n, b, k, fp[k - 1 + offset] + 1, fp[k + 1 + offset]);
        } // for
        fp[delta + offset] = onp_snake(m, a, n, b, delta, fp[delta - 1 + offset] + 1, fp[delta + 1 + offset]);
        p += 1;
    } // while

    allocator.allocate(allocator.context, fp, size * sizeof(*fp), 0);
    return delta + 2 * (p - 1);
} // onp_compare


/*
static inline size_t
onp_snake_matches(size_t const m, char const a[static restrict m],
    size_t const n, char const b[static restrict n],
    intmax_t const k,
    intmax_t const lower,
    intmax_t const upper,
    size_t const p,
    GVA_Matches matches)
{
    size_t col = MAX(lower, upper);
    size_t row = col - k;

    intmax_t const d_row = m - row;
    intmax_t const d_col = n - col;

    while (row < m && col < n && a[row] == b[col])
    {
        row += 1;
        col += 1;
        size_t const lcs_pos = (row + col - (n - m) - 2 * p + ABS(d_row - d_col)) / 2;
        matches.uniq[lcs_pos - 1] = 0;
        if (lcs_pos > matches.max_lcs_pos)
        {
            matches.max_lcs_pos = lcs_pos;
            matches.uniq[lcs_pos - 1] = 1;
            matches.matches[lcs_pos - 1].row = row - 1;
            matches.matches[lcs_pos - 1].col = col - 1;
        } // if
    } // while
    return col;
} // onp_snake_matches


static GVA_Matches
onp_compare_matches(GVA_Allocator const allocator,
    size_t const m, char const a[static restrict m],
    size_t const n, char const b[static restrict n])
{
    intmax_t const delta = n - m;
    size_t const offset = m + 1;
    size_t const size = m + n + 3;
    intmax_t* const restrict fp = allocator.allocate(allocator.context, NULL, 0, size * sizeof(*fp));
    if (fp == NULL)
    {
        return (GVA_Matches) {NULL};
    } // if

    for (intmax_t i = -1; i <= delta; ++i)
    {
        fp[i + offset] = -1;
    } // for

    GVA_Matches result =
    {
        .matches = allocator.allocate(allocator.context, NULL, 0, sizeof(*result.matches) * m),
        .uniq = allocator.allocate(allocator.context, NULL, 0, sizeof(*result.uniq) * m),
        .max_lcs_pos = 0,
    };
    if (result.matches == NULL || result.uniq == NULL)
    {
        result.uniq = allocator.allocate(allocator.context, result.uniq, sizeof(*result.uniq) * m, 0);
        result.matches = allocator.allocate(allocator.context, result.matches, sizeof(*result.matches) * m, 0);
        allocator.allocate(allocator.context, fp, size * sizeof(*fp), 0);
        return (GVA_Matches) {NULL};
    } // if

    size_t p = 0;
    while ((size_t) fp[delta + offset] != n)
    {
        fp[-p - 1 + offset] = -1;
        fp[delta + p + 1 + offset] = -1;
        for (intmax_t k = -p; k <= delta - 1; ++k)
        {
            fp[k + offset] = onp_snake_matches(m, a, n, b, k, fp[k - 1 + offset] + 1, fp[k + 1 + offset], p, result);
        } // for
        for (intmax_t k = delta + p; k >= delta + 1; --k)
        {
            fp[k + offset] = onp_snake_matches(m, a, n, b, k, fp[k - 1 + offset] + 1, fp[k + 1 + offset], p, result);
        } // for
        fp[delta + offset] = onp_snake_matches(m, a, n, b, delta, fp[delta - 1 + offset] + 1, fp[delta + 1 + offset], p, result);
        p += 1;
    } // while

    allocator.allocate(allocator.context, fp, size * sizeof(*fp), 0);
    result.distance = delta + 2 * (p - 1);
    return result;
} // onp_compare_matches
*/


inline size_t
gva_edit_distance(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs])
{
    return len_ref > len_obs ?
           onp_compare(allocator, len_obs, observed, len_ref, reference) :
           onp_compare(allocator, len_ref, reference, len_obs, observed);
} // gva_edit_distance


/*
inline GVA_Matches
gva_edit_distance_matches(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs])
{
    if (len_ref > len_obs)
    {
        GVA_Matches result = onp_compare_matches(allocator, len_obs, observed, len_ref, reference);
        for (size_t i = 0; i < result.max_lcs_pos; ++i)
        {
            gva_uint const temp = result.matches[i].row;
            result.matches[i].row = result.matches[i].col;
            result.matches[i].col = temp;
        } // for
        return result;
    } // if

    return onp_compare_matches(allocator, len_ref, reference, len_obs, observed);
} // gva_edit_distance_matches
*/
