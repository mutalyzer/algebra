#include <stdint.h>     // uint8_t
#include <stddef.h>     // NULL, size_t
#include <string.h>     // memcpy, strncmp


#include "../include/allocator.h"   // GVA_Allocator
#include "../include/edit.h"        // gva_edit_distance
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_lcs_graph_*, gva_edges
#include "../include/relations.h"   // GVA_Relation, GVA_CONTAINS, GVA_DISJOINT,
                                    // GVA_EQUIVALENT GVA_IS_CONTAINED, GVA_OVERLAP
#include "../include/variant.h"     // GVA_Variant
#include "array.h"      // ARRAY_DESTROY, array_length
#include "common.h"     // GVA_NULL, gva_uint, MAX, MIN
#include "bitset.h"     // bitset_*


static uint8_t const NUC_A = 0x1;
static uint8_t const NUC_C = 0x2;
static uint8_t const NUC_G = 0x4;
static uint8_t const NUC_T = 0x8;


static inline uint8_t
nucleotides(size_t const len, char const sequence[static len])
{
    static uint8_t const MASK[256] =
    {
        ['A'] = NUC_A,
        ['C'] = NUC_C,
        ['G'] = NUC_G,
        ['T'] = NUC_T,
    };
    static uint8_t const UNIVERSE = 0xF;

    uint8_t mask = 0x0;
    for (size_t i = 0; mask < UNIVERSE && i < len; ++i)
    {
        mask |= MASK[(size_t) sequence[i]];
    } // for
    return mask;
} // nucleotides


static size_t
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
        if (graph.nodes[i].row > end)
        {
            continue;
        } // if

        for (gva_uint j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
        {
            if (graph.nodes[graph.edges[j].tail].row + graph.nodes[graph.edges[j].tail].length < start)
            {
                continue;
            } // if

            GVA_Variant variant;
            gva_uint const count = gva_edges(graph.observed.str,
                                             graph.nodes[i], graph.nodes[graph.edges[j].tail],
                                             i == graph.source, graph.nodes[graph.edges[j].tail].edges == GVA_NULL,
                                             &variant);

            if (variant.end > variant.start)
            {
                bitset_add(dels, variant.start - offset, variant.end + count - 1 - offset);
            } // if

            uint8_t const mask = nucleotides(variant.sequence.len, variant.sequence.str);
            if ((mask & NUC_A) == NUC_A)
            {
                bitset_add(as, variant.start - offset, variant.end + count - offset);
            } // if
            if ((mask & NUC_C) == NUC_C)
            {
                bitset_add(cs, variant.start - offset, variant.end + count - offset);
            } // if
            if ((mask & NUC_G) == NUC_G)
            {
                bitset_add(gs, variant.start - offset, variant.end + count - offset);
            } // if
            if ((mask & NUC_T) == NUC_T)
            {
                bitset_add(ts, variant.start - offset, variant.end + count - offset);
            } // if
        } // for
    } // for
    return counter;
} // bitset_fill


size_t disjoint_count = 0;


GVA_Relation
gva_compare(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref],
    GVA_Variant const lhs, GVA_Variant const rhs)
{
    GVA_Relation relation;

    if (lhs.start == rhs.start && lhs.end == rhs.end &&
        lhs.sequence.len == rhs.sequence.len &&
        strncmp(lhs.sequence.str, rhs.sequence.str, MIN(lhs.sequence.len, rhs.sequence.len)) == 0)
    {
        return GVA_EQUIVALENT;
    } // if

    if (lhs.start > rhs.end || rhs.start > lhs.end)
    {
        return GVA_DISJOINT;
    } // if

    size_t const start = MIN(lhs.start, rhs.start);
    size_t const end = MAX(lhs.end, rhs.end);

    size_t const lhs_len = (lhs.start - start) + lhs.sequence.len + (end - lhs.end);
    char* lhs_obs = allocator.allocate(allocator.context, NULL, 0, lhs_len);
    if (lhs_len > 0 && lhs_obs == NULL)
    {
        return -1;
    } // if

    size_t const rhs_len = (rhs.start - start) + rhs.sequence.len + (end - rhs.end);
    char* rhs_obs = allocator.allocate(allocator.context, NULL, 0, rhs_len);
    if (rhs_len > 0 && rhs_obs == NULL)
    {
        lhs_obs = allocator.allocate(allocator.context, lhs_obs, lhs_len, 0);
        return -1;
    } // if

    size_t const lhs_dist = gva_edit_distance(allocator, lhs.end - lhs.start, reference + lhs.start, lhs.sequence.len, lhs.sequence.str);
    size_t const rhs_dist = gva_edit_distance(allocator, rhs.end - rhs.start, reference + rhs.start, rhs.sequence.len, rhs.sequence.str);

    memcpy(lhs_obs, reference + start, lhs.start - start);
    memcpy(lhs_obs + lhs.start - start, lhs.sequence.str, lhs.sequence.len);
    memcpy(lhs_obs + lhs.start - start + lhs.sequence.len, reference + lhs.end, end - lhs.end);

    memcpy(rhs_obs, reference + start, rhs.start - start);
    memcpy(rhs_obs + rhs.start - start, rhs.sequence.str, rhs.sequence.len);
    memcpy(rhs_obs + rhs.start - start + rhs.sequence.len, reference + rhs.end, end - rhs.end);

    size_t const dist = gva_edit_distance(allocator, lhs_len, lhs_obs, rhs_len, rhs_obs);

    if (lhs_dist + rhs_dist == dist)
    {
        relation = GVA_DISJOINT;
    } // if
    else if (lhs_dist - rhs_dist == dist)
    {
        relation = GVA_CONTAINS;
    } // if
    else if (rhs_dist - lhs_dist == dist)
    {
        relation = GVA_IS_CONTAINED;
    } // if
    else
    {
        GVA_LCS_Graph lhs_graph = gva_lcs_graph_init(allocator, end - start, reference + start, lhs_len, lhs_obs, start);
        GVA_LCS_Graph rhs_graph = gva_lcs_graph_init(allocator, end - start, reference + start, rhs_len, rhs_obs, start);

        gva_uint const len = end - start + 1;
        size_t* lhs_dels = bitset_init(allocator, len); // can be one shorter
        size_t* lhs_as = bitset_init(allocator, len);
        size_t* lhs_cs = bitset_init(allocator, len);
        size_t* lhs_gs = bitset_init(allocator, len);
        size_t* lhs_ts = bitset_init(allocator, len);

        size_t* rhs_dels = bitset_init(allocator, len); // can be one shorter
        size_t* rhs_as = bitset_init(allocator, len);
        size_t* rhs_cs = bitset_init(allocator, len);
        size_t* rhs_gs = bitset_init(allocator, len);
        size_t* rhs_ts = bitset_init(allocator, len);

        size_t const start_intersection = MAX(lhs.start, rhs.start);
        size_t const end_intersection = MIN(lhs.end, rhs.end);

        // could be done on intersection instead of union
        bitset_fill(lhs_graph, start, start_intersection, end_intersection, lhs_dels, lhs_as, lhs_cs, lhs_gs, lhs_ts);
        gva_lcs_graph_destroy(allocator, lhs_graph);

        bitset_fill(rhs_graph, start, start_intersection, end_intersection, rhs_dels, rhs_as, rhs_cs, rhs_gs, rhs_ts);
        gva_lcs_graph_destroy(allocator, rhs_graph);

        size_t const common =
            bitset_intersection_cnt(lhs_dels, rhs_dels) +
            bitset_intersection_cnt(lhs_as, rhs_as) +
            bitset_intersection_cnt(lhs_cs, rhs_cs) +
            bitset_intersection_cnt(lhs_gs, rhs_gs) +
            bitset_intersection_cnt(lhs_ts, rhs_ts);

        if (common > 0)
        {
            relation = GVA_OVERLAP;
        } // if
        else
        {
            relation = GVA_DISJOINT;
            disjoint_count += 1;
        } // else

        lhs_ts = bitset_destroy(allocator, lhs_ts);
        lhs_gs = bitset_destroy(allocator, lhs_gs);
        lhs_cs = bitset_destroy(allocator, lhs_cs);
        lhs_as = bitset_destroy(allocator, lhs_as);
        lhs_dels = bitset_destroy(allocator, lhs_dels);

        rhs_ts = bitset_destroy(allocator, rhs_ts);
        rhs_gs = bitset_destroy(allocator, rhs_gs);
        rhs_cs = bitset_destroy(allocator, rhs_cs);
        rhs_as = bitset_destroy(allocator, rhs_as);
        rhs_dels = bitset_destroy(allocator, rhs_dels);
    } // else

    rhs_obs = allocator.allocate(allocator.context, rhs_obs, rhs_len, 0);
    lhs_obs = allocator.allocate(allocator.context, lhs_obs, lhs_len, 0);

    return relation;
} // gva_compare
