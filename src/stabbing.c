#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t


#include <stdio.h>      // DEBUG


#include "../include/allocator.h"   // GVA_Allocator
#include "../include/stabbing.h"    // GVA_Stabbing_Entry, GVA_Stabbing_Index, gva_stabbing_index_*
#include "array.h"      // ARRAY_*, array_*
#include "common.h"     // GVA_NULL, gva_uint, MAX


GVA_Stabbing_Index
gva_stabbing_index_init(GVA_Allocator const allocator, size_t const len_ref, GVA_Stabbing_Entry* const entries)
{
    GVA_Stabbing_Index index = {entries, NULL, len_ref};

    index.start_table = allocator.allocate(allocator.context, NULL, 0, index.len_ref * sizeof(*index.start_table));
    gva_uint* stack = allocator.allocate(allocator.context, NULL, 0, index.len_ref * sizeof(*stack));
    if (index.start_table == NULL || stack == NULL)
    {
        stack = allocator.allocate(allocator.context, NULL, index.len_ref * sizeof(*stack), 0);
        index.entries = ARRAY_DESTROY(allocator, index.entries);
        index.start_table = allocator.allocate(allocator.context, index.start_table, index.len_ref * sizeof(*index.start_table), 0);
        index.len_ref = 0;
        return index;
    } // if

    gva_uint prev = GVA_NULL;
    gva_uint idx = 0;
    gva_uint top = GVA_NULL;
    for (gva_uint i = 0; i < index.len_ref; ++i)
    {
        index.start_table[i] = GVA_NULL;
        //fprintf(stderr, "index iteration: %u\n", i);

        while (top != GVA_NULL && entries[stack[top]].end < i)
        {
            //fprintf(stderr, "Popping %u [%u--%u] from stack\n", top, entries[stack[top]].start, entries[stack[top]].end);
            top -= 1;
        } // while

        while (idx < array_length(entries) &&
               i >= entries[idx].start && i <= entries[idx].end &&
               !(top != GVA_NULL && entries[idx].start == entries[stack[top]].start && entries[idx].end == entries[stack[top]].end))
        {
            do
            {
                //fprintf(stderr, "Pushing %u [%u--%u] on stack\n", idx, entries[idx].start, entries[idx].end);
                top += 1;
                stack[top] = idx;

                entries[idx].parent = prev;
                while (entries[idx].parent != GVA_NULL && entries[entries[idx].parent].end < entries[idx].end)
                {
                    entries[idx].parent = entries[entries[idx].parent].parent;
                } // while

                if (entries[idx].parent != GVA_NULL)
                {
                    entries[idx].leftsibling = entries[entries[idx].parent].rightchild;
                    entries[entries[idx].parent].rightchild = idx;
                } // if
                else if (idx != stack[0])
                {
                    entries[idx].leftsibling = stack[0];
                } // else
                prev = idx;

                idx += 1;
            } while (idx < array_length(entries) &&
                     entries[idx].start == entries[stack[top]].start && entries[idx].end == entries[stack[top]].end);
        } // while

        if (top != GVA_NULL)
        {
            //fprintf(stderr, "start_table[%u] := %d\n", i, top);
            index.start_table[i] = stack[top];
        } // if
    } // for
    stack = allocator.allocate(allocator.context, stack, index.len_ref * sizeof(*stack), 0);

    return index;
} // gva_stabbing_index_init


void
gva_stabbing_index_destroy(GVA_Allocator const allocator, GVA_Stabbing_Index index)
{
    index.entries = ARRAY_DESTROY(allocator, index.entries);
    index.start_table = allocator.allocate(allocator.context, index.start_table, index.len_ref * sizeof(*index.start_table), 0);
    index.len_ref = 0;
} // gva_stabbing_index_destroy


gva_uint*
gva_stabbing_index_stab(GVA_Allocator const allocator, GVA_Stabbing_Index const index,
    gva_uint const start, gva_uint const end)
{
    gva_uint* results = NULL;
    gva_uint* stack = NULL;
    gva_uint stab = index.start_table[end];
    while (stab != GVA_NULL)
    {
        //fprintf(stderr, "stab: %2u: [%2u, %2u]\n", stab, index.entries[stab].start, index.entries[stab].end);
        ARRAY_APPEND(allocator, results, stab);

        if (index.entries[stab].leftsibling != GVA_NULL && start <= index.entries[index.entries[stab].leftsibling].end)
        {
            ARRAY_APPEND(allocator, stack, index.entries[stab].leftsibling);
            while (array_length(stack) > 0)
            {
                array_header(stack)->length -= 1;
                gva_uint trav = stack[array_length(stack)];

                //fprintf(stderr, "trav: %2u: [%2u, %2u]\n", trav, index.entries[trav].start, index.entries[trav].end);
                ARRAY_APPEND(allocator, results, trav);

                if (index.entries[trav].leftsibling != GVA_NULL && start <= index.entries[index.entries[trav].leftsibling].end)
                {
                    ARRAY_APPEND(allocator, stack, index.entries[trav].leftsibling);
                } // if

                if (index.entries[trav].rightchild != GVA_NULL && start <= index.entries[index.entries[trav].rightchild].end)
                {
                    ARRAY_APPEND(allocator, stack, index.entries[trav].rightchild);
                } // if
            } // while
        } // if
        stab = index.entries[stab].parent;
    } // while
    stack = ARRAY_DESTROY(allocator, stack);

    return results;
} // gva_stabbing_index_stab
