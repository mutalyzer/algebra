#include <inttypes.h>   // intmax_t
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <string.h>     // memset

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "align.h"      // LCS_Alignment, LCS_Node, lcs_align
#include "array.h"      // ARRAY_APPEND
#include "common.h"     // ABS, MAX, MIN


typedef struct
{
    GVA_Allocator const allocator;
    size_t const len_ref;
    char const* const restrict reference;
    size_t const len_obs;
    char const* const restrict observed;
    size_t const shift;
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
            gva_uint const len = ARRAY_APPEND(context.allocator, lcs->nodes, ((LCS_Node) {
                match_row + context.shift, match_col, length, false, 0, GVA_NULL, GVA_NULL
            })) - 1;
            if (lcs->index[lcs_pos].head != GVA_NULL)
            {
                lcs->nodes[lcs->index[lcs_pos].tail].next = len;
            } // if
            else
            {
                lcs->index[lcs_pos].head = len;
            } // else
            lcs->index[lcs_pos].tail = len;

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
        gva_uint const len = ARRAY_APPEND(context.allocator, lcs->nodes, ((LCS_Node) {
            match_row + context.shift, match_col, length, false, 0, GVA_NULL, GVA_NULL
        })) - 1;
        if (lcs->index[lcs_pos].head != GVA_NULL)
        {
            lcs->nodes[lcs->index[lcs_pos].tail].next = len;
        } // if
        else
        {
            lcs->index[lcs_pos].head = len;
        } // else
        lcs->index[lcs_pos].tail = len;

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
    size_t const shift)
{
    intmax_t const delta = len_obs - len_ref;
    size_t const offset = len_ref + 1;
    size_t const size = len_ref + len_obs + 3;

    LCS_Alignment lcs = {0, NULL, NULL};

    lcs.index = allocator.allocate(allocator.context, NULL, 0, MIN(len_ref, len_obs) * sizeof(*lcs.index));
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
        shift,
        .diagonals = allocator.allocate(allocator.context, NULL, 0, size * sizeof(*context.diagonals)),
    };
    if (context.diagonals == NULL)
    {
        lcs.index = allocator.allocate(allocator.context, lcs.index, MIN(len_ref, len_obs) * sizeof(*lcs.index), lcs.length * sizeof(*lcs.index));
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
    while (context.diagonals[offset + delta] <= len)
    {
        for (intmax_t idx = lower - p; idx < delta; ++idx)
        {
            context.diagonals[offset + idx] = expand(context, idx, p, &lcs);
        } // for
        for (intmax_t idx = upper + p; idx > delta; --idx)
        {
            context.diagonals[offset + idx] = expand(context, idx, p, &lcs);
        } // for
        context.diagonals[offset + delta] = expand(context, delta, p, &lcs);

        p += 1;
    } // while

    lcs.index = allocator.allocate(allocator.context, lcs.index, MIN(len_ref, len_obs) * sizeof(*lcs.index), lcs.length * sizeof(*lcs.index));
    allocator.allocate(allocator.context, context.diagonals, size * sizeof(*context.diagonals), 0);
    return lcs;
} // lcs_align
