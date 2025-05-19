#include <inttypes.h>   // imaxabs
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, ptrdiff_t, size_t

#include "../include/alloc.h"   // VA_Allocator
#include "../include/array.h"   // va_array_*
#include "../include/edit.h"


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
    VA_Allocator const allocator;
    size_t const len_ref;
    char const* const restrict reference;
    size_t const len_obs;
    char const* const restrict observed;
    size_t* const restrict diagonals;
    uint32_t* const restrict lcs_index;
} Context;


static size_t
expand(Context const context,
       ptrdiff_t const idx,
       size_t const p,
       size_t* const len_lcs,
       VA_LCS_Node2** lcs_nodes)
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
            if (lcs_pos + 1 > *len_lcs)
            {
                *len_lcs = lcs_pos + 1;
            } // if
            size_t const length = row - match_row;
            va_array_append(context.allocator, *lcs_nodes, ((VA_LCS_Node2) {match_row, match_col, length, lcs_pos, context.lcs_index[lcs_pos], -1}));
            context.lcs_index[lcs_pos] = va_array_length(*lcs_nodes) - 1;
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
        if (lcs_pos + 1 > *len_lcs)
        {
            *len_lcs = lcs_pos + 1;
        } // if
        size_t const length = row - match_row;
        va_array_append(context.allocator, *lcs_nodes, ((VA_LCS_Node2) {match_row, match_col, length, lcs_pos, context.lcs_index[lcs_pos], -1}));
        context.lcs_index[lcs_pos] = va_array_length(*lcs_nodes) - 1;
    } // if

    return steps;
} // expand


size_t
va_edit2(VA_Allocator const allocator,
         size_t const len_ref,
         char const reference[static restrict len_ref],
         size_t const len_obs,
         char const observed[static restrict len_obs],
         VA_LCS_Node2** restrict lcs_nodes,
         uint32_t** restrict lcs_index)
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
        .diagonals = allocator.alloc(allocator.context, NULL, 0, size * sizeof(*context.diagonals)),
        .lcs_index = allocator.alloc(allocator.context, NULL, 0, umin(len_ref, len_obs) * sizeof(*context.lcs_index)),
    };
    if (context.diagonals == NULL || context.lcs_index == NULL)
    {
        allocator.alloc(allocator.context, context.diagonals, size * sizeof(*context.diagonals), 0);
        allocator.alloc(allocator.context, context.lcs_index, umin(len_ref, len_obs) * sizeof(*context.lcs_index), 0);
        return -1;
    } // if

    for (size_t i = 0; i < size; ++i)
    {
        context.diagonals[i] = 0;
    } // for

    for (size_t i = 0; i < umin(len_ref, len_obs); ++i)
    {
        context.lcs_index[i] = -1;
    } // for

    *lcs_nodes = NULL;

    size_t const lower = delta > 0 ? 0 : delta;
    size_t const upper = delta > 0 ? delta : 0;

    size_t const len = umax(len_ref, len_obs) - imaxabs(delta);
    size_t p = 0;
    size_t max_lcs_pos = 0;
    while (context.diagonals[offset + delta] <= len)
    {
        for (ptrdiff_t idx = lower - p; idx < delta; ++idx)
        {
            context.diagonals[offset + idx] = expand(context, idx, p, &max_lcs_pos, lcs_nodes);
        } // for
        for (ptrdiff_t idx = upper + p; idx > delta; --idx)
        {
            context.diagonals[offset + idx] = expand(context, idx, p, &max_lcs_pos, lcs_nodes);
        } // for
        context.diagonals[offset + delta] = expand(context, delta, p, &max_lcs_pos, lcs_nodes);

        p += 1;
    } // while

    allocator.alloc(allocator.context, context.diagonals, size * sizeof(*context.diagonals), 0);

    // transfer ownership
    *lcs_index = allocator.alloc(allocator.context, context.lcs_index, umin(len_ref, len_obs) * sizeof(*context.lcs_index), max_lcs_pos * sizeof(*context.lcs_index));

    return max_lcs_pos;
} // va_edit2
