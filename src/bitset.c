// FIXME: __builtin_popcountll
#include <limits.h>     // CHAR_BIT
#include <stddef.h>     // NULL, size_t
#include <string.h>     // memset


#include "../include/allocator.h"   // GVA_Allocator
#include "array.h"      // ARRAY_*, array_*
#include "bitset.h"     // bitset_*


inline size_t*
bitset_init(GVA_Allocator const allocator, size_t const size)
{
    size_t const len = (size + sizeof(size_t) * CHAR_BIT - 1) / (sizeof(size_t) * CHAR_BIT);
    size_t* bitset = array_init(allocator, len, sizeof(*bitset));
    if (bitset == NULL)
    {
        return NULL;
    } // if
    array_header(bitset)->length = len;
    return memset(bitset, 0, len * sizeof(*bitset));
} // bitset_init


inline void
bitset_add(size_t bitset[static 1], size_t const value)
{
    bitset[value / (sizeof(size_t) * CHAR_BIT)] |= (1ULL << (value % (sizeof(size_t) * CHAR_BIT)));
} // bitset_add


inline size_t
bitset_intersection_cnt(size_t const lhs[static 1], size_t const rhs[static 1])
{
    size_t const len = array_length(lhs);
    size_t count = 0;
    for (size_t i = 0; i < len; ++i)
    {
        count += __builtin_popcountll(lhs[i] & rhs[i]);
    } // for
    return count;
} // bitset_intersection_cnt


inline size_t*
bitset_destroy(GVA_Allocator const allocator, size_t bitset[static 1])
{
    return ARRAY_DESTROY(allocator, bitset);
} // bitset_destroy


size_t
bitset_fill(GVA_LCS_Graph const graph,
    gva_uint const offset,
    gva_uint const start, gva_uint const end,
    size_t dels[static restrict 1],
    size_t as[static restrict 1],
    size_t cs[static restrict 1],
    size_t gs[static restrict 1],
    size_t ts[static restrict 1])
{
    size_t counter = 0;
    for (size_t i = 0; i < array_length(graph.nodes); ++i)
    {
        for (gva_uint j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
        {
            if (graph.nodes[i].row > end || graph.nodes[graph.edges[j].tail].row + graph.nodes[graph.edges[j].tail].length < start)
            {
                continue;
            } // if

            GVA_Variant variant;
            gva_uint const count = gva_edges(graph.observed.str,
                                             graph.nodes[i], graph.nodes[graph.edges[j].tail],
                                             i == graph.source, graph.nodes[graph.edges[j].tail].edges == GVA_NULL,
                                             &variant);

            for (gva_uint z = 0; z < count; ++z)
            {
                for (gva_uint k = variant.start; k < variant.end; ++k)
                {
                    bitset_add(dels, k + z - offset);
                    counter += 1;
                } // for
            } // for

            if (memchr(variant.sequence.str, 'A', variant.sequence.len))
            {
                for (gva_uint k = variant.start; k < variant.end + count; ++k)
                {
                    bitset_add(as, k - offset);
                    counter += 1;
                } // for
            } // if
            if (memchr(variant.sequence.str, 'C', variant.sequence.len))
            {
                for (gva_uint k = variant.start; k < variant.end + count; ++k)
                {
                    bitset_add(cs, k - offset);
                    counter += 1;
                } // for
            } // if
            if (memchr(variant.sequence.str, 'G', variant.sequence.len))
            {
                for (gva_uint k = variant.start; k < variant.end + count; ++k)
                {
                    bitset_add(gs, k - offset);
                    counter += 1;
                } // for
            } // if
            if (memchr(variant.sequence.str, 'T', variant.sequence.len))
            {
                for (gva_uint k = variant.start; k < variant.end + count; ++k)
                {
                    bitset_add(ts, k - offset);
                    counter += 1;
                } // for
            } // if
        } // for
    } // for
    return counter;
} // bitfill
