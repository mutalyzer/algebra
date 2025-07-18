#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <string.h>     // memset


#include <stdio.h>      // DEBUG


#include "../include/allocator.h"   // GVA_Allocator
#include "../include/stabbing.h"    // GVA_Stabbing_Entry, GVA_Stabbing_Index, gva_stabbing_index_*
#include "array.h"      // ARRAY_*, array_*
#include "common.h"     // GVA_NULL, gva_uint, MAX


GVA_Stabbing_Index
gva_stabbing_index_init(GVA_Allocator const allocator, size_t const len_ref)
{
    gva_uint* start_table = allocator.allocate(allocator.context, NULL, 0, len_ref * sizeof(*start_table));  // OVERFLOW
    if (start_table == NULL)
    {
        return (GVA_Stabbing_Index) {0, NULL, NULL};
    } // if

    memset(start_table, 0, len_ref * sizeof(*start_table));
    GVA_Stabbing_Entry* entries = NULL;
    ARRAY_APPEND(allocator, entries, ((GVA_Stabbing_Entry) {GVA_NULL, GVA_NULL, GVA_NULL, 0, len_ref}));

    return (GVA_Stabbing_Index) {len_ref, start_table, entries};
} // gva_stabbing_index_init


inline size_t
gva_stabbing_index_add(GVA_Allocator const allocator, GVA_Stabbing_Index index[static 1],
    size_t const start, size_t const end)
{
    return ARRAY_APPEND(allocator, index->entries, ((GVA_Stabbing_Entry) {GVA_NULL, GVA_NULL, GVA_NULL, start, end}));
} // gva_stabbing_index_add


void
gva_stabbing_index_build(GVA_Allocator const allocator, GVA_Stabbing_Index const index)
{
    gva_uint const len = array_length(index.entries);
    if (len <= 1)
    {
        return;
    } // if

    gva_uint* stack = allocator.allocate(allocator.context, NULL, 0, index.len_ref * sizeof(*stack));
    if (stack == NULL)
    {
        return;
    } // if

    gva_uint prev = 0;
    gva_uint idx = 1;
    gva_uint top = 0;
    stack[top] = 0;
    for (gva_uint i = 0; i < index.len_ref; ++i)
    {
        fprintf(stderr, "index iteration: %u\n", i);

        while (index.entries[stack[top]].end < i)
        {
            fprintf(stderr, "Popping [%u, %u] from stack\n", index.entries[stack[top]].start, index.entries[stack[top]].end);
            top -= 1;
        } // while

        while (idx < len && i >= index.entries[idx].start && i <= index.entries[idx].end &&
            !(index.entries[idx].start == index.entries[stack[top]].start && index.entries[idx].end == index.entries[stack[top]].end))
        {
            do
            {
                fprintf(stderr, "Pushing [%u, %u] on stack\n", index.entries[idx].start, index.entries[idx].end);
                top += 1;
                stack[top] = idx;

                index.entries[idx].parent = prev;
                while (index.entries[idx].parent != GVA_NULL && index.entries[index.entries[idx].parent].end < index.entries[idx].end)
                {
                    index.entries[idx].parent = index.entries[index.entries[idx].parent].parent;
                } // while
                index.entries[idx].leftsibling = index.entries[index.entries[idx].parent].rightchild;
                index.entries[index.entries[idx].parent].rightchild = idx;
                prev = idx;

                idx += 1;
           } while (idx < len && index.entries[idx].start == index.entries[stack[top]].start && index.entries[idx].end == index.entries[stack[top]].end);
        } /// while

        fprintf(stderr, "start_table[%u] := [%u, %u]\n", i, index.entries[stack[top]].start, index.entries[stack[top]].end);
        index.start_table[i] = stack[top];
    } // for

    stack = allocator.allocate(allocator.context, stack, index.len_ref * sizeof(*stack), 0);
} // gva_stabbing_index_init


void
gva_stabbing_index_destroy(GVA_Allocator const allocator, GVA_Stabbing_Index index[static 1])
{
    index->entries = ARRAY_DESTROY(allocator, index->entries);
    index->start_table = allocator.allocate(allocator.context, index->start_table, index->len_ref * sizeof(*index->start_table), 0);
    index->len_ref = 0;
} // gva_stabbing_index_destroy


gva_uint*
gva_stabbing_index_intersect(GVA_Allocator const allocator, GVA_Stabbing_Index const index,
    gva_uint const start, gva_uint const end)
{
    gva_uint* results = NULL;
    gva_uint* stack = NULL;
    gva_uint stab = index.start_table[end];
    while (stab != GVA_NULL && index.entries[stab].parent != GVA_NULL)
    {
        fprintf(stderr, "stab: %2u: [%2u, %2u]\n", stab, index.entries[stab].start, index.entries[stab].end);
        ARRAY_APPEND(allocator, results, stab);

        if (index.entries[stab].leftsibling != GVA_NULL && start <= index.entries[index.entries[stab].leftsibling].end)
        {
            ARRAY_APPEND(allocator, stack, index.entries[stab].leftsibling);
            while (array_length(stack) > 0)
            {
                array_header(stack)->length -= 1;
                gva_uint trav = stack[array_length(stack)];

                fprintf(stderr, "trav: %2u: [%2u, %2u]\n", trav, index.entries[trav].start, index.entries[trav].end);
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

    if (index.entries[stab].rightchild != GVA_NULL && start <= index.entries[index.entries[stab].rightchild].end)
    {
        ARRAY_APPEND(allocator, stack, index.entries[stab].rightchild);
        while (array_length(stack) > 0)
        {
            array_header(stack)->length -= 1;
            gva_uint trav = stack[array_length(stack)];

            fprintf(stderr, "Rtrav: %2u: [%2u, %2u]\n", trav, index.entries[trav].start, index.entries[trav].end);
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


    stack = ARRAY_DESTROY(allocator, stack);

    return results;
} // gva_stabbing_index_intersect
