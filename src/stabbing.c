#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <string.h>     // memset

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/stabbing.h"    // GVA_Stabbing_Entry, GVA_Stabbing_Index, gva_stabbing_index_*
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "array.h"      // ARRAY_*, array_*
#include "common.h"     // MAX


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
    ARRAY_APPEND(allocator, entries, ((GVA_Stabbing_Entry) {GVA_NULL, GVA_NULL, GVA_NULL, 0, len_ref, 0, 0, 0}));

    return (GVA_Stabbing_Index) {len_ref, start_table, entries};
} // gva_stabbing_index_init


void
gva_stabbing_index_destroy(GVA_Allocator const allocator, GVA_Stabbing_Index index[static 1])
{
    index->entries = ARRAY_DESTROY(allocator, index->entries);
    index->start_table = allocator.allocate(allocator.context, index->start_table, index->len_ref * sizeof(*index->start_table), 0);
    index->len_ref = 0;
} // gva_stabbing_index_destroy


inline size_t
gva_stabbing_index_add(GVA_Allocator const allocator, GVA_Stabbing_Index index[static 1],
    gva_uint const start, gva_uint const end, gva_uint const inserted,
    gva_uint const distance, size_t const data)
{
    return ARRAY_APPEND(allocator, index->entries, ((GVA_Stabbing_Entry) {GVA_NULL, GVA_NULL, GVA_NULL, start, end, inserted, distance, data}));
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
        while (index.entries[stack[top]].end < i)
        {
            top -= 1;
        } // while

        while (idx < len && i >= index.entries[idx].start && i <= index.entries[idx].end &&
            !(index.entries[idx].start == index.entries[stack[top]].start && index.entries[idx].end == index.entries[stack[top]].end))
        {
            do
            {
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

        index.start_table[i] = stack[top];
    } // for

    stack = allocator.allocate(allocator.context, stack, index.len_ref * sizeof(*stack), 0);
} // gva_stabbing_index_init


gva_uint*
gva_stabbing_index_self_intersect(GVA_Allocator const allocator, GVA_Stabbing_Index const index)
{
    gva_uint* results = NULL;
    for (size_t i = 1; i < array_length(index.entries); ++i)
    {
        for (size_t j = i + 1; j < array_length(index.entries); ++j)
        {
            if (index.entries[i].end < index.entries[j].start)
            {
                break;
            } // if
            if (index.entries[j].end >= index.entries[i].start && index.entries[j].start <= index.entries[i].end)
            {
                ARRAY_APPEND(allocator, results, 0);
            } // if
        } // for
    } // for
    return results;
} // gva_stabbing_index_self_intersect


gva_uint*
gva_stabbing_index_intersect(GVA_Allocator const allocator, GVA_Stabbing_Index const index,
    gva_uint const start, gva_uint const end)
{
    gva_uint* results = NULL;
    gva_uint* stack = NULL;
    gva_uint stab = index.start_table[end];

    if (index.entries[stab].rightchild != GVA_NULL)
    {
        ARRAY_APPEND(allocator, stack, index.entries[stab].rightchild);
        while (array_length(stack) > 0)
        {
            array_header(stack)->length -= 1;
            gva_uint trav = stack[array_length(stack)];
            if (end >= index.entries[trav].start)
            {
                ARRAY_APPEND(allocator, results, trav);
            } // if

            if (index.entries[trav].leftsibling != GVA_NULL)
            {
                ARRAY_APPEND(allocator, stack, index.entries[trav].leftsibling);
            } // if

            if (index.entries[trav].rightchild != GVA_NULL && end >= index.entries[index.entries[trav].rightchild].start)
            {
                ARRAY_APPEND(allocator, stack, index.entries[trav].rightchild);
            } // if
        } // while
    } // if

    while (stab != GVA_NULL)
    {
        if (stab > 0)
        {
            ARRAY_APPEND(allocator, results, stab);
        } // if

        if (index.entries[stab].leftsibling != GVA_NULL && start <= index.entries[index.entries[stab].leftsibling].end)
        {
            ARRAY_APPEND(allocator, stack, index.entries[stab].leftsibling);
            while (array_length(stack) > 0)
            {
                array_header(stack)->length -= 1;
                gva_uint trav = stack[array_length(stack)];
                ARRAY_APPEND(allocator, results, trav);

                if (index.entries[trav].leftsibling != GVA_NULL && start <= index.entries[index.entries[trav].leftsibling].end)
                {
                    ARRAY_APPEND(allocator, stack, index.entries[trav].leftsibling);
                } // if

                if (index.entries[trav].rightchild != GVA_NULL)
                {
                    ARRAY_APPEND(allocator, stack, index.entries[trav].rightchild);
                } // if
            } // while
        } // if

        stab = index.entries[stab].parent;
    } // while
    stack = ARRAY_DESTROY(allocator, stack);

    return results;
} // gva_stabbing_index_intersect
