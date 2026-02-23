#include <stdbool.h>    // bool
#include <stddef.h>     // NULL, size_t
#include <string.h>     // memcpy

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/compare.h"     // gva_compare_graphs
#include "../include/edit.h"        // gva_edit_distance
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_lcs_graph_*, gva_edges
#include "../include/relations.h"   // GVA_Relation, GVA_CONTAINS, GVA_DISJOINT,
                                    // GVA_EQUIVALENT GVA_IS_CONTAINED, GVA_OVERLAP
#include "../include/variant.h"     // GVA_Variant, gva_variant_eq
#include "common.h"     // MAX, MIN
#include "bitset.h"     // bitset_*


GVA_Relation
gva_compare_graphs(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref],
    GVA_LCS_Graph const lhs, GVA_LCS_Graph const rhs)
{
    if (gva_variant_eq(lhs.supremal, rhs.supremal))
    {
        return GVA_EQUIVALENT;
    } // if

    if (lhs.supremal.start > rhs.supremal.end || rhs.supremal.start > lhs.supremal.end)
    {
        return GVA_DISJOINT;
    } // if

    size_t const start = MIN(lhs.supremal.start, rhs.supremal.start);
    size_t const end = MAX(lhs.supremal.end, rhs.supremal.end);

    size_t const len_lhs = (lhs.supremal.start - start) + lhs.supremal.sequence.len + (end - lhs.supremal.end);
    size_t const len_rhs = (rhs.supremal.start - start) + rhs.supremal.sequence.len + (end - rhs.supremal.end);

    size_t distance = 0;
    if (len_lhs == 0)
    {
        distance = len_rhs;
    } // if
    else if (len_rhs == 0)
    {
        distance = len_lhs;
    } // if
    else
    {
        GVA_String observed_lhs = gva_string_init(allocator, len_lhs);
        GVA_String observed_rhs = gva_string_init(allocator, len_rhs);
        if (observed_lhs.str == NULL || observed_rhs.str == NULL)
        {
            gva_string_destroy(allocator, observed_rhs);
            gva_string_destroy(allocator, observed_lhs);
            return GVA_DISJOINT;  // FIXME: OOM
        } // if

        memcpy((char*) observed_lhs.str, reference + start, lhs.supremal.start - start);
        memcpy((char*) observed_lhs.str + lhs.supremal.start - start, lhs.supremal.sequence.str, lhs.supremal.sequence.len);
        memcpy((char*) observed_lhs.str + lhs.supremal.start - start + lhs.supremal.sequence.len, reference + lhs.supremal.end, end - lhs.supremal.end);

        memcpy((char*) observed_rhs.str, reference + start, rhs.supremal.start - start);
        memcpy((char*) observed_rhs.str + rhs.supremal.start - start, rhs.supremal.sequence.str, rhs.supremal.sequence.len);
        memcpy((char*) observed_rhs.str + rhs.supremal.start - start + rhs.supremal.sequence.len, reference + rhs.supremal.end, end - rhs.supremal.end);

        distance = gva_edit_distance(allocator, observed_lhs.len, observed_lhs.str, observed_rhs.len, observed_rhs.str);
        gva_string_destroy(allocator, observed_rhs);
        gva_string_destroy(allocator, observed_lhs);
    } // else

    if (gva_lcs_graph_distance(lhs) + gva_lcs_graph_distance(rhs) == distance)
    {
        return GVA_DISJOINT;
    } // if

    if (gva_lcs_graph_distance(lhs) - gva_lcs_graph_distance(rhs) == distance)
    {
        return GVA_CONTAINS;
    } // if

    if (gva_lcs_graph_distance(rhs) - gva_lcs_graph_distance(lhs) == distance)
    {
        return GVA_IS_CONTAINED;
    } // if

    size_t const len = end - start + 1;
    size_t* lhs_dels = bitset_init(allocator, len);
    size_t* lhs_as = bitset_init(allocator, len);
    size_t* lhs_cs = bitset_init(allocator, len);
    size_t* lhs_gs = bitset_init(allocator, len);
    size_t* lhs_ts = bitset_init(allocator, len);

    size_t* rhs_dels = bitset_init(allocator, len);
    size_t* rhs_as = bitset_init(allocator, len);
    size_t* rhs_cs = bitset_init(allocator, len);
    size_t* rhs_gs = bitset_init(allocator, len);
    size_t* rhs_ts = bitset_init(allocator, len);

    size_t const start_intersection = MAX(lhs.supremal.start, rhs.supremal.start);
    size_t const end_intersection = MIN(lhs.supremal.end, rhs.supremal.end);

    gva_lcs_graph_uniq_atomics(lhs, start, start_intersection, end_intersection, lhs_dels, lhs_as, lhs_cs, lhs_gs, lhs_ts);
    gva_lcs_graph_uniq_atomics(rhs, start, start_intersection, end_intersection, rhs_dels, rhs_as, rhs_cs, rhs_gs, rhs_ts);

    bool const overlap = bitset_intersection_cnt(lhs_dels, rhs_dels) > 0 ||
        bitset_intersection_cnt(lhs_as, rhs_as) > 0 ||
        bitset_intersection_cnt(lhs_cs, rhs_cs) > 0 ||
        bitset_intersection_cnt(lhs_gs, rhs_gs) > 0 ||
        bitset_intersection_cnt(lhs_ts, rhs_ts) > 0;

    rhs_ts = bitset_destroy(allocator, rhs_ts);
    rhs_gs = bitset_destroy(allocator, rhs_gs);
    rhs_cs = bitset_destroy(allocator, rhs_cs);
    rhs_as = bitset_destroy(allocator, rhs_as);
    rhs_dels = bitset_destroy(allocator, rhs_dels);

    lhs_ts = bitset_destroy(allocator, lhs_ts);
    lhs_gs = bitset_destroy(allocator, lhs_gs);
    lhs_cs = bitset_destroy(allocator, lhs_cs);
    lhs_as = bitset_destroy(allocator, lhs_as);
    lhs_dels = bitset_destroy(allocator, lhs_dels);

    return overlap ? GVA_OVERLAP : GVA_DISJOINT;
} // gva_compare_graphs
