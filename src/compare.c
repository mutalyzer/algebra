#include <stdio.h>
#include <string.h>                 // memchr, memcpy, strncmp
#include "../include/allocator.h"   // GVA_Allocator
#include "../include/edit.h"        // gva_edit_distance
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_lcs_graph_*, gva_edges
#include "../include/relations.h"   // GVA_Relation, GVA_DISJOINT
#include "../include/variant.h"     // GVA_Variant
#include "array.h"      // ARRAY_DESTROY, array_length
#include "common.h"     // GVA_NULL, gva_uint, MAX, MIN
#include "bitset.h"     // bitset_*



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
    }

    if (lhs.start > rhs.end || rhs.start > lhs.end)
    {
        return GVA_DISJOINT;
    }

    size_t const start = MIN(lhs.start, rhs.start);
    size_t const end = MAX(lhs.end, rhs.end);
//    fprintf(stderr, "start: %zu end: %zu\n", start, end);

    size_t const lhs_len = (lhs.start - start) + lhs.sequence.len + (end - lhs.end);
//    fprintf(stderr, "lhs_len: %zu\n", lhs_len);
    char* lhs_obs = allocator.allocate(allocator.context, NULL, 0, lhs_len);
    if (lhs_len > 0 && lhs_obs == NULL)
    {
        return -1;
    }

    size_t const rhs_len = (rhs.start - start) + rhs.sequence.len + (end - rhs.end);
//    fprintf(stderr, "rhs_len: %zu\n", rhs_len);
    char* rhs_obs = allocator.allocate(allocator.context, NULL, 0, rhs_len);
    if (rhs_len > 0 && rhs_obs == NULL)
    {
        lhs_obs = allocator.allocate(allocator.context, lhs_obs, lhs_len, 0);
        return -1;
    }

    size_t const lhs_dist = gva_edit_distance(allocator, lhs.end - lhs.start, reference + lhs.start, lhs.sequence.len, lhs.sequence.str);
    size_t const rhs_dist = gva_edit_distance(allocator, rhs.end - rhs.start, reference + rhs.start, rhs.sequence.len, rhs.sequence.str);

    memcpy(lhs_obs, reference + start, lhs.start - start);
    memcpy(lhs_obs + lhs.start - start, lhs.sequence.str, lhs.sequence.len);
    memcpy(lhs_obs + lhs.start - start + lhs.sequence.len, reference + lhs.end, end - lhs.end);

//    fprintf(stderr, "lhs_obs(%zu): %.*s\n", lhs_len, (int) lhs_len, lhs_obs);

    memcpy(rhs_obs, reference + start, rhs.start - start);
    memcpy(rhs_obs + rhs.start - start, rhs.sequence.str, rhs.sequence.len);
    memcpy(rhs_obs + rhs.start - start + rhs.sequence.len, reference + rhs.end, end - rhs.end);

//    fprintf(stderr, "rhs_obs(%zu): %.*s\n", rhs_len, (int) rhs_len, rhs_obs);

    size_t const dist = gva_edit_distance(allocator, lhs_len, lhs_obs, rhs_len, rhs_obs);

//    fprintf(stderr, "distances: %zu, %zu, %zu\n", lhs_dist, rhs_dist, dist);

    if (lhs_dist + rhs_dist == dist)
    {
        relation = GVA_DISJOINT;
    } //
    else if (lhs_dist - rhs_dist == dist)
    {
        relation = GVA_CONTAINS;
    } // if
    else if (rhs_dist - lhs_dist == dist)
    {
        relation = GVA_IS_CONTAINED;
    } // if
    else {

        // inject random bits
        GVA_LCS_Graph lhs_graph = gva_lcs_graph_init(allocator, end - start, reference + start, lhs_len, lhs_obs, start);
        GVA_LCS_Graph rhs_graph = gva_lcs_graph_init(allocator, end - start, reference + start, rhs_len, rhs_obs, start);

        gva_uint lhs_sup_len = end - start;
//        fprintf(stderr, "lhs suplen %d\n", lhs_sup_len);
        size_t* lhs_dels = bitset_init(allocator, lhs_sup_len + 1); // can be one shorter
        size_t* lhs_as = bitset_init(allocator, lhs_sup_len + 1);
        size_t* lhs_cs = bitset_init(allocator, lhs_sup_len + 1);
        size_t* lhs_gs = bitset_init(allocator, lhs_sup_len + 1);
        size_t* lhs_ts = bitset_init(allocator, lhs_sup_len + 1);

        gva_uint rhs_sup_len = end - start;
//        fprintf(stderr, "rhs suplen %d\n", rhs_sup_len);
        size_t* rhs_dels = bitset_init(allocator, rhs_sup_len + 1); // can be one shorter
        size_t* rhs_as = bitset_init(allocator, rhs_sup_len + 1);
        size_t* rhs_cs = bitset_init(allocator, rhs_sup_len + 1);
        size_t* rhs_gs = bitset_init(allocator, rhs_sup_len + 1);
        size_t* rhs_ts = bitset_init(allocator, rhs_sup_len + 1);

        size_t const start_accent = MAX(lhs.start, rhs.start);
        size_t const end_accent = MIN(lhs.end, rhs.end);

        // could be done on intersection instead of union
        bitset_fill(lhs_graph, start, start_accent, end_accent, lhs_dels, lhs_as, lhs_cs, lhs_gs, lhs_ts);
        gva_lcs_graph_destroy(allocator, lhs_graph);

//        for (size_t i = 0; i < array_length(lhs_dels); ++i)
//        {
//            fprintf(stderr, "%2zu: %016zx\n", i, lhs_dels[i]);
//            fprintf(stderr, "%2zu: %016zx\n", i, lhs_as[i]);
//            fprintf(stderr, "%2zu: %016zx\n", i, lhs_cs[i]);
//            fprintf(stderr, "%2zu: %016zx\n", i, lhs_gs[i]);
//            fprintf(stderr, "%2zu: %016zx\n", i, lhs_ts[i]);
//        } // for

        bitset_fill(rhs_graph, start, start_accent, end_accent, rhs_dels, rhs_as, rhs_cs, rhs_gs, rhs_ts);
        gva_lcs_graph_destroy(allocator, rhs_graph);

//        for (size_t i = 0; i < array_length(rhs_dels); ++i)
//        {
//            fprintf(stderr, "%2zu: %016zx\n", i, rhs_dels[i]);
//            fprintf(stderr, "%2zu: %016zx\n", i, rhs_as[i]);
//            fprintf(stderr, "%2zu: %016zx\n", i, rhs_cs[i]);
//            fprintf(stderr, "%2zu: %016zx\n", i, rhs_gs[i]);
//            fprintf(stderr, "%2zu: %016zx\n", i, rhs_ts[i]);
//        } // for

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
}