#include <inttypes.h>   // intmax_t
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <string.h>     // memset

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "align.h"      // LCS_Alignment, LCS_Matches, LCS_Node, lcs_align*
#include "array.h"      // ARRAY_APPEND
#include "common.h"     // ABS, MAX, MIN


#include <stdio.h>  // DEBUG


typedef struct
{
    GVA_Allocator const allocator;
    size_t const len_ref;
    char const* const restrict reference;
    size_t const len_obs;
    char const* const restrict observed;
    size_t const offset;
    size_t* const restrict diagonals;
} Context;


static size_t
expand(Context const context,
    intmax_t const idx,
    size_t const p,
    LCS_Alignment lcs[static 1])
{
    intmax_t const delta = context.len_obs - context.len_ref;
    size_t const offset = context.len_ref + 1;

    size_t row;
    size_t col;
    size_t end;
    if (idx > 0)
    {
        row = context.diagonals[offset + idx];
        col = row + idx;
        end = MAX(context.diagonals[offset + idx - 1] - 1, context.diagonals[offset + idx + 1]);
    } // if
    else if (idx < 0)
    {
        col = context.diagonals[offset + idx];
        row = col - idx;
        end = MAX(context.diagonals[offset + idx - 1], context.diagonals[offset + idx + 1] - 1);
    } // if
    else
    {
        row = context.diagonals[offset + idx];
        col = row + idx;
        end = MAX(context.diagonals[offset + idx - 1], context.diagonals[offset + idx + 1]);
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
            intmax_t const d_row = context.len_ref - row;
            intmax_t const d_col = context.len_obs - col;
            size_t const lcs_pos = (row + col - ABS(delta) - 2 * p + ABS(d_row - d_col)) / 2 - 1;
            size_t const length = row - match_row;
            gva_uint const idx = ARRAY_APPEND(context.allocator, lcs->nodes,
                ((LCS_Node)
                {
                    .match = {match_row + context.offset, match_col, length},
                    .idx = GVA_NULL,
                    .next = GVA_NULL,
                }));
            if (lcs->index[lcs_pos].head != GVA_NULL)
            {
                lcs->nodes[lcs->index[lcs_pos].tail].next = idx;
            } // if
            else
            {
                lcs->index[lcs_pos].head = idx;
            } // else
            lcs->index[lcs_pos].tail = idx;

            if (lcs_pos + 1 > lcs->length)
            {
                lcs->length = lcs_pos + 1;
            } // if
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
        intmax_t const d_row = context.len_ref - row;
        intmax_t const d_col = context.len_obs - col;
        size_t const lcs_pos = (row + col - ABS(delta) - 2 * p + ABS(d_row - d_col)) / 2 - 1;
        size_t const length = row - match_row;
        gva_uint const idx = ARRAY_APPEND(context.allocator, lcs->nodes,
            ((LCS_Node)
            {
                .match = {match_row + context.offset, match_col, length},
                .idx = GVA_NULL,
                .next = GVA_NULL,
            }));
        if (lcs->index[lcs_pos].head != GVA_NULL)
        {
            lcs->nodes[lcs->index[lcs_pos].tail].next = idx;
        } // if
        else
        {
            lcs->index[lcs_pos].head = idx;
        } // else
        lcs->index[lcs_pos].tail = idx;

        if (lcs_pos + 1 > lcs->length)
        {
            lcs->length = lcs_pos + 1;
        } // if
    } // if
    return steps;
} // expand


LCS_Alignment
lcs_align(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs],
    size_t const offset)
{
    intmax_t const delta = len_obs - len_ref;
    size_t const size = len_ref + len_obs + 3;

    LCS_Alignment lcs = {0};

    lcs.index = allocator.allocate(allocator.context, NULL,
        0, MIN(len_ref, len_obs) * sizeof(*lcs.index));
    if (lcs.index == NULL)
    {
        return lcs;
    } // if

    Context const context =
    {
        allocator,
        len_ref,
        reference,
        len_obs,
        observed,
        offset,
        .diagonals = allocator.allocate(allocator.context, NULL,
            0, size * sizeof(*context.diagonals)),
    };
    if (context.diagonals == NULL)
    {
        lcs.index = allocator.allocate(allocator.context, lcs.index,
            MIN(len_ref, len_obs) * sizeof(*lcs.index), lcs.length * sizeof(*lcs.index));
        return lcs;
    } // if

    memset(context.diagonals, 0, size * sizeof(*context.diagonals));

    for (size_t i = 0; i < MIN(len_ref, len_obs); ++i)
    {
        lcs.index[i].head = GVA_NULL;
        lcs.index[i].tail = GVA_NULL;
        lcs.index[i].count = 0;
    } // for

    size_t const lower = delta > 0 ? 0 : delta;
    size_t const upper = delta > 0 ? delta : 0;
    size_t const len = MAX(len_ref, len_obs) - ABS(delta);
    size_t p = 0;
    while (context.diagonals[len_ref + 1 + delta] <= len)
    {
        for (intmax_t idx = lower - p; idx < delta; ++idx)
        {
            context.diagonals[len_ref + 1 + idx] = expand(context, idx, p, &lcs);
        } // for
        for (intmax_t idx = upper + p; idx > delta; --idx)
        {
            context.diagonals[len_ref + 1 + idx] = expand(context, idx, p, &lcs);
        } // for
        context.diagonals[len_ref + 1 + delta] = expand(context, delta, p, &lcs);

        p += 1;
    } // while

    lcs.index = allocator.allocate(allocator.context, lcs.index,
        MIN(len_ref, len_obs) * sizeof(*lcs.index), lcs.length * sizeof(*lcs.index));
    allocator.allocate(allocator.context, context.diagonals,
        size * sizeof(*context.diagonals), 0);
    return lcs;
} // lcs_align


static inline size_t
onp_snake(size_t const m, char const a[static restrict m],
    size_t const n, char const b[static restrict n],
    intmax_t const k, intmax_t const lower, intmax_t const upper,
    intmax_t const delta, size_t const p, bool const swapped,
    LCS_Matches result[static restrict 1])
{
    size_t col = MAX(lower, upper);
    size_t row = col - k;

    if (swapped)
    {
        fprintf(stderr, "  snake start: (%zu, %zu)\n", col, row);
    } // if
    else
    {
        fprintf(stderr, "  snake start: (%zu, %zu)\n", row, col);
    } // else

    intmax_t const d_row = m - row;
    intmax_t const d_col = n - col;

    while (row < m && col < n && a[row] == b[col])
    {
        row += 1;
        col += 1;
        size_t const lcs_pos = (row + col - ABS(delta) - 2 * p + ABS(d_row - d_col)) / 2;

        if (swapped)
        {
            fprintf(stderr, "    snake match @ %zu: (%zu, %zu)\n", lcs_pos - 1, col - 1, row - 1);
        } // if
        else
        {
            fprintf(stderr, "    snake match @ %zu: (%zu, %zu)\n", lcs_pos - 1, row - 1, col - 1);
        } // else

        result->uniq[lcs_pos - 1] = 0;
        if (lcs_pos > result->max_lcs_pos)
        {
            result->max_lcs_pos = lcs_pos;
            result->uniq[lcs_pos - 1] = 1;
            if (swapped)
            {
                result->match[lcs_pos - 1].row = col - 1;
                result->match[lcs_pos - 1].col = row - 1;
            } // if
            else
            {
                result->match[lcs_pos - 1].row = row - 1;
                result->match[lcs_pos - 1].col = col - 1;
            } // else
        } // if
    } // while

    if (swapped)
    {
        fprintf(stderr, "  snake end: (%zu, %zu)\n", col, row);
    } // if
    else
    {
        fprintf(stderr, "  snake end: (%zu, %zu)\n", row, col);
    } // else

    return col;
} // onp_snake


static LCS_Matches
onp_compare(GVA_Allocator const allocator,
    size_t const m, char const a[static restrict m],
    size_t const n, char const b[static restrict n],
    bool const swapped)
{
    intmax_t const delta = n - m;
    size_t const offset = m + 1;
    size_t const size = m + n + 3;

    intmax_t* restrict fp = allocator.allocate(allocator.context, NULL, 0, size * sizeof(*fp));
    LCS_Matches result =
    {
        .match = allocator.allocate(allocator.context, NULL, 0, sizeof(*result.match) * m),
        .uniq = allocator.allocate(allocator.context, NULL, 0, sizeof(*result.uniq) * m),
    };
    if (fp == NULL || result.match == NULL || result.uniq == NULL)
    {
        result.uniq = allocator.allocate(allocator.context, result.uniq, sizeof(*result.uniq) * m, 0);
        result.match = allocator.allocate(allocator.context, result.match, sizeof(*result.match) * m, 0);
        fp = allocator.allocate(allocator.context, fp, size * sizeof(*fp), 0);
        return (LCS_Matches) {NULL};
    } // if

    for (intmax_t i = -1; i <= delta; ++i)
    {
        fp[i + offset] = -1;
    } // for

    size_t p = 0;
    while ((size_t) fp[delta + offset] != n)
    {
        fprintf(stderr, "p: %zu\n", p);
        fp[-p - 1 + offset] = -1;
        fp[delta + p + 1 + offset] = -1;
        for (intmax_t k = -p; k <= delta - 1; ++k)
        {
            fp[k + offset] = onp_snake(m, a, n, b, k, fp[k - 1 + offset] + 1, fp[k + 1 + offset], delta, p, swapped, &result);
        } // for
        for (intmax_t k = delta + p; k >= delta + 1; --k)
        {
            fp[k + offset] = onp_snake(m, a, n, b, k, fp[k - 1 + offset] + 1, fp[k + 1 + offset], delta, p, swapped, &result);
        } // for
        fp[delta + offset] = onp_snake(m, a, n, b, delta, fp[delta - 1 + offset] + 1, fp[delta + 1 + offset], delta, p, swapped, &result);
        p += 1;
    } // while

    allocator.allocate(allocator.context, fp, size * sizeof(*fp), 0);

    result.distance = delta + 2 * (p - 1);
    return result;
} // onp_compare


inline LCS_Matches
lcs_align_one(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs])
{
    return len_ref > len_obs ?
           onp_compare(allocator, len_obs, observed, len_ref, reference, true) :
           onp_compare(allocator, len_ref, reference, len_obs, observed, false);
} // lcs_align_one
