#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // intmax_t
#include <string.h>     // memcpy

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/edit.h"        // gva_edit_distance
#include "../include/index.h"       // GVA_Index, GVA_Query_Result, gva_index_*
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_lcs_graph_*
#include "../include/relations.h"   // GVA_Relation, GVA_CONTAINS,
                                    // GVA_DISJOINT, GVA_EQUIVALENT,
                                    // GVA_IS_CONTAINED
#include "../include/string.h"      // GVA_String
#include "../include/types.h"       // GVA_NULL, GVA_Interval, gva_uint
#include "../include/variant.h"     // GVA_Variant, gva_variant_dup
#include "array.h"              // ARRAY_*, array_length
#include "bitset.h"             // bitset_*
#include "common.h"             // ABS, MAX, MIN
#include "hash_table.h"         // NOT_FOUND, HASH_TABLE_*, hash_table_*
#include "interval_tree.h"      // Interval_Tree, interval_tree_*
#include "trie.h"               // Trie, trie_*


#include <stdio.h>      // DEBUG


typedef struct
{
    gva_uint link;  // X-OR link between interval_idx and allele_idx
    gva_uint next;
} Join;


typedef struct
{
    gva_uint     id_idx;
    gva_uint     distance;
    GVA_Interval join;
} Allele;


struct GVA_Index
{
    GVA_Allocator allocator;
    GVA_String    reference;
    Interval_Tree intervals;
    Trie          inserted;
    Trie          ids;
    Allele*       alleles;
    Join*         join;
};


static inline GVA_Relation
relation_from_included(size_t const included,
    size_t const lhs_distance, size_t const rhs_distance,
    gva_uint excluded[static 1])
{
    *excluded = lhs_distance + rhs_distance - 2 * included;
    if (*excluded == 0)
    {
        return GVA_EQUIVALENT;
    } // if
    if (lhs_distance == included)
    {
        return GVA_IS_CONTAINED;
    } // if
    if (rhs_distance == included)
    {
         return GVA_CONTAINS;
    } // if
    return GVA_OVERLAP;
} // relation_from_included


static inline GVA_Variant
variant_from_index(GVA_Index const self[static 1], size_t const idx)
{
    return (GVA_Variant)
    {
        self->intervals.nodes[idx].start,
        self->intervals.nodes[idx].end,
        trie_string(self->inserted, self->intervals.nodes[idx].inserted)
    };
} // variant_from_index


static size_t
variants_distance(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref],
    GVA_Variant const lhs, GVA_Variant const rhs)
{
    size_t const start = MIN(lhs.start, rhs.start);
    size_t const end = MAX(lhs.end, rhs.end);

    size_t const len_lhs = (lhs.start - start) + lhs.sequence.len + (end - lhs.end);
    size_t const len_rhs = (rhs.start - start) + rhs.sequence.len + (end - rhs.end);

    if (len_lhs == 0)
    {
        return len_rhs;
    } // if
    if (len_rhs == 0)
    {
        return len_lhs;
    } // if

    GVA_String observed_lhs = gva_string_init(allocator, len_lhs);
    GVA_String observed_rhs = gva_string_init(allocator, len_rhs);
    if (observed_lhs.str == NULL || observed_rhs.str == NULL)
    {
        gva_string_destroy(allocator, observed_rhs);
        gva_string_destroy(allocator, observed_lhs);
        return -1;  // FIXME: OOM
    } // if

    memcpy((char*) observed_lhs.str, reference + start, lhs.start - start);
    if (lhs.sequence.str != NULL)
    {
        memcpy((char*) observed_lhs.str + lhs.start - start, lhs.sequence.str, lhs.sequence.len);
    } // if
    memcpy((char*) observed_lhs.str + lhs.start - start + lhs.sequence.len, reference + lhs.end, end - lhs.end);

    memcpy((char*) observed_rhs.str, reference + start, rhs.start - start);
    if (rhs.sequence.str != NULL)
    {
        memcpy((char*) observed_rhs.str + rhs.start - start, rhs.sequence.str, rhs.sequence.len);
    } // if
    memcpy((char*) observed_rhs.str + rhs.start - start + rhs.sequence.len, reference + rhs.end, end - rhs.end);

    size_t const distance = gva_edit_distance(allocator, observed_lhs.len, observed_lhs.str, observed_rhs.len, observed_rhs.str);

    gva_string_destroy(allocator, observed_rhs);
    gva_string_destroy(allocator, observed_lhs);

    return distance;
} // variants_distance


static size_t
variants_included(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref],
    GVA_Variant const lhs, GVA_Variant const rhs,
    size_t* const excluded)
{
    size_t const start = MIN(lhs.start, rhs.start);
    size_t const end = MAX(lhs.end, rhs.end);

    GVA_String observed_lhs = gva_string_init(allocator, (lhs.start - start) + lhs.sequence.len + (end - lhs.end));
    GVA_String observed_rhs = gva_string_init(allocator, (rhs.start - start) + rhs.sequence.len + (end - rhs.end));

    memcpy((char*) observed_lhs.str, reference + start, lhs.start - start);
    if (lhs.sequence.str != NULL)
    {
        memcpy((char*) observed_lhs.str + lhs.start - start, lhs.sequence.str, lhs.sequence.len);
    } // if
    memcpy((char*) observed_lhs.str + lhs.start - start + lhs.sequence.len, reference + lhs.end, end - lhs.end);

    memcpy((char*) observed_rhs.str, reference + start, rhs.start - start);
    if (rhs.sequence.str != NULL)
    {
        memcpy((char*) observed_rhs.str + rhs.start - start, rhs.sequence.str, rhs.sequence.len);
    } // if
    memcpy((char*) observed_rhs.str + rhs.start - start + rhs.sequence.len, reference + rhs.end, end - rhs.end);

    GVA_LCS_Graph lhs_graph = gva_lcs_graph_init(allocator, end - start, reference + start, observed_lhs.len, observed_lhs.str, start);
    GVA_LCS_Graph rhs_graph = gva_lcs_graph_init(allocator, end - start, reference + start, observed_rhs.len, observed_rhs.str, start);

    gva_uint const len = end - start + 1;
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

    size_t const start_intersection = MAX(lhs.start, rhs.start);
    size_t const end_intersection = MIN(lhs.end, rhs.end);

    gva_lcs_graph_uniq_atomics(lhs_graph, start, start_intersection, end_intersection, lhs_dels, lhs_as, lhs_cs, lhs_gs, lhs_ts);
    gva_lcs_graph_destroy(allocator, lhs_graph, false);

    gva_lcs_graph_uniq_atomics(rhs_graph, start, start_intersection, end_intersection, rhs_dels, rhs_as, rhs_cs, rhs_gs, rhs_ts);
    gva_lcs_graph_destroy(allocator, rhs_graph, false);

    size_t const included =
        bitset_intersection_cnt(lhs_dels, rhs_dels) +
        bitset_intersection_cnt(lhs_as, rhs_as) +
        bitset_intersection_cnt(lhs_cs, rhs_cs) +
        bitset_intersection_cnt(lhs_gs, rhs_gs) +
        bitset_intersection_cnt(lhs_ts, rhs_ts);

    *excluded =
        bitset_intersection_cnt(lhs_dels, lhs_dels) +
        bitset_intersection_cnt(rhs_dels, rhs_dels) +
        bitset_intersection_cnt(lhs_as, lhs_as) +
        bitset_intersection_cnt(rhs_as, rhs_as) +
        bitset_intersection_cnt(lhs_cs, lhs_cs) +
        bitset_intersection_cnt(rhs_cs, rhs_cs) +
        bitset_intersection_cnt(lhs_gs, lhs_gs) +
        bitset_intersection_cnt(rhs_gs, rhs_gs) +
        bitset_intersection_cnt(lhs_ts, lhs_ts) +
        bitset_intersection_cnt(rhs_ts, rhs_ts);

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

    gva_string_destroy(allocator, observed_rhs);
    gva_string_destroy(allocator, observed_lhs);

    return included;
} // variants_included


static inline size_t
variants_with_distance(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref],
    GVA_Variant const lhs, size_t const distance_lhs,
    GVA_Variant const rhs, size_t const distance_rhs)
{
    size_t const distance = variants_distance(allocator, len_ref, reference, lhs, rhs);

    if (distance == 0 || distance_rhs - distance_lhs == distance)
    {
        return distance_lhs;
    } // if
    if (distance_lhs - distance_rhs == distance)
    {
        return distance_rhs;
    } // if

    size_t excluded = 0;
    size_t const included = variants_included(allocator, len_ref, reference, lhs, rhs, &excluded);
    if (included == 0)
    {
        return 0;  // disjoint
    } // if
    return MAX(1, (double) (2 * included - 1) / excluded * MIN(distance_lhs, distance_rhs));
} // variants_with_distance


inline GVA_Index*
gva_index_init(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref])
{
    GVA_Index* index = allocator.allocate(allocator.context, NULL, 0, sizeof(*index));
    if (index == NULL)
    {
        return NULL;
    } // if

    index->allocator = allocator;
    index->reference = (GVA_String) {len_ref, reference};

    index->intervals = interval_tree_init(allocator);
    index->inserted = trie_init(allocator);

    index->ids = trie_init(allocator);

    index->alleles = NULL;
    index->join = NULL;

    return index;
} // gva_index_init


inline GVA_Index*
gva_index_destroy(GVA_Index* const self)
{
    if (self == NULL)
    {
        return NULL;
    } // if

    interval_tree_destroy(&self->intervals);
    trie_destroy(&self->inserted);
    trie_destroy(&self->ids);
    self->alleles = ARRAY_DESTROY(self->allocator, self->alleles);
    self->join = ARRAY_DESTROY(self->allocator, self->join);

    return self->allocator.allocate(self->allocator.context, self, sizeof(*self), 0);
} // gva_index_destroy


void
gva_index_insert(GVA_Index* restrict const self,
    size_t const len_id, char const id[static restrict len_id],
    GVA_Variant const variant, size_t const distance)
{
    size_t const id_idx = trie_insert(&self->ids, len_id, id);
    size_t allele_idx = array_length(self->alleles) - 1;

    // new allele
    if (array_length(self->alleles) == 0 ||
        self->alleles[allele_idx].id_idx != id_idx)
    {
        allele_idx = ARRAY_APPEND(self->allocator, self->alleles,
            ((Allele)
            {
                .id_idx = id_idx,
                .join.start = array_length(self->join),
            }));
    } // if

    // add variant
    gva_uint const inserted_idx = trie_insert(&self->inserted,
        variant.sequence.len, variant.sequence.str);
    gva_uint const tmp_idx = ARRAY_APPEND(self->allocator, self->intervals.nodes,
        ((Interval_Tree_Node)
        {
            .child = {GVA_NULL, GVA_NULL},
            .start = variant.start,
            .end = variant.end,
            .max = variant.end,
            .inserted = inserted_idx,
            .alleles = GVA_NULL,
            .distance = distance,
        }));
    gva_uint const node_idx = interval_tree_insert(&self->intervals, tmp_idx);
    // undo append: interval already in the tree
    if (node_idx != tmp_idx)
    {
        array_header(self->intervals.nodes)->length -= 1;
    } // if
    self->intervals.nodes[node_idx].alleles = ARRAY_APPEND(self->allocator, self->join,
        ((Join)
        {
            .link = node_idx ^ allele_idx,
            .next = self->intervals.nodes[node_idx].alleles,
        }));

    // update allele
    self->alleles[allele_idx].distance += distance;
    self->alleles[allele_idx].join.end = self->intervals.nodes[node_idx].alleles + 1;
} // gva_index_insert


GVA_Query_Result
gva_index_query(GVA_Allocator const allocator,
    GVA_Index* const self, GVA_LCS_Graph const graph)
{
    static size_t const INITIAL_SIZE = 1024;

    struct Allele_Entry
    {
        HASH_TABLE_KEY;     // index into self->alleles
        gva_uint included;
        gva_uint head;      // singly linked list in hits
        gva_uint tail;
    }* entries = hash_table_init(allocator, INITIAL_SIZE, sizeof(*entries));

    struct Hit
    {
        GVA_Interval join;
        GVA_Interval query;
        gva_uint     included;  // double purpose: distance or included
        gva_uint     excluded;
        GVA_Relation relation;
        gva_uint     next;
    }* hits = NULL;

    // Phase 1: for every local supremal part:
    //          - find candidates based on interval query
    //          - calculate the distance with the candidates,
    //            and discard disjoint candidates early based on distance
    for (size_t i = 0; i < array_length(graph.dom_nodes) - 1; ++i)
    {
        GVA_Variant const variant = gva_lcs_graph_local_supremal(graph, i, i + 1);

        gva_uint* intervals = interval_tree_intersection(allocator, self->intervals, variant.start, variant.end);
        for (size_t j = 0; j < array_length(intervals); ++j)
        {
            gva_uint const node_idx = intervals[j];

            size_t const distance = variants_distance(allocator, self->reference.len, self->reference.str,
                variant_from_index(self, node_idx), variant);
            if (self->intervals.nodes[node_idx].distance + graph.dom_nodes[i + 1].distance - graph.dom_nodes[i].distance <= distance)
            {
                // Disjoint based on distance
                continue;
            } // if

            for (gva_uint k = self->intervals.nodes[node_idx].alleles; k != GVA_NULL; k = self->join[k].next)
            {
                gva_uint const tail = ARRAY_APPEND(allocator, hits,
                    ((struct Hit)
                    {
                        .join = {k, k + 1},
                        .query = {i, i + 1},
                        .included = distance,
                        .next = GVA_NULL,
                    }));

                gva_uint const allele_idx = self->join[k].link ^ node_idx;
                size_t const idx = HASH_TABLE_INDEX(entries, allele_idx);
                if (entries[idx].gva_key != allele_idx)
                {
                    HASH_TABLE_SET(allocator, entries, allele_idx,
                        ((struct Allele_Entry)
                        {
                            .gva_key = allele_idx,
                            .head = tail,
                            .tail = tail,
                        }));
                    continue;
                } // if

                hits[entries[idx].tail].next = tail;
                entries[idx].tail = tail;
            } // for
        } // for
        intervals = ARRAY_DESTROY(allocator, intervals);
    } // for

    gva_uint* ranges = NULL;

    // Phase 2
    for (size_t idx = 0; idx < array_header(entries)->capacity; ++idx)
    {
        if (entries[idx].gva_key == NOT_FOUND)
        {
            continue;
        } // if

        //fprintf(stderr, GVA_STRING_FMT "\n", GVA_STRING_PRINT(gva_index_id(self, entries[idx].gva_key)));
        gva_uint prev = GVA_NULL;
        for (gva_uint i = entries[idx].head; i != GVA_NULL; i = hits[i].next)
        {
            //fprintf(stderr, "  %2u %2u %2u\n", hits[i].join.start, hits[i].query.start, hits[i].included);

            // look aheads for multi hits: build new local supremal
            size_t start_ranges = array_length(ranges);
            GVA_Variant supremal = {0};
            size_t distance = 0;

            // check single join hit with multiple query hits
            gva_uint start = hits[i].query.start;
            gva_uint end = hits[i].query.end;
            while (hits[i].next != GVA_NULL && hits[i].join.start == hits[hits[i].next].join.start)
            {
                if (start_ranges == array_length(ranges))
                {
                    ARRAY_APPEND(allocator, ranges, hits[i].query.start);
                } // if
                i = hits[i].next;
                ARRAY_APPEND(allocator, ranges, hits[i].query.start);
                if (hits[i].query.start > end)
                {
                    GVA_Variant const variant = gva_lcs_graph_local_supremal(graph, start, end);
                    if (distance == 0)
                    {
                        supremal = gva_variant_dup(allocator, variant);
                    } // if
                    else
                    {
                        supremal.sequence = gva_string_concat(allocator, supremal.sequence,
                            (GVA_String) {variant.start - supremal.end, self->reference.str + supremal.end});
                        supremal.sequence = gva_string_concat(allocator, supremal.sequence, variant.sequence);
                        supremal.end = variant.end;
                    } // else
                    distance += graph.dom_nodes[end].distance - graph.dom_nodes[start].distance;
                    start = hits[i].query.start;
                } // if
                end = hits[i].query.end;
                //fprintf(stderr, "   \" %2u %2u (query)\n", hits[i].query.start, hits[i].included);
            } // while
            if (start_ranges < array_length(ranges))
            {
                GVA_Variant const variant = gva_lcs_graph_local_supremal(graph, start, end);
                if (distance == 0)
                {
                    supremal = gva_variant_dup(allocator, variant);
                } // if
                else
                {
                    supremal.sequence = gva_string_concat(allocator, supremal.sequence,
                        (GVA_String) {variant.start - supremal.end, self->reference.str + supremal.end});
                    supremal.sequence = gva_string_concat(allocator, supremal.sequence, variant.sequence);
                    supremal.end = variant.end;
                } // else
                distance += (graph.dom_nodes[end].distance - graph.dom_nodes[start].distance);
                //fprintf(stderr, "  SOLVE multiple query: " GVA_VARIANT_FMT_SPDI " (%zu)\n", GVA_VARIANT_PRINT_SPDI("NC_000001.11", supremal), distance);
                gva_uint const node_idx = self->join[hits[i].join.start].link ^ entries[idx].gva_key;
                hits[i].included = variants_with_distance(allocator, self->reference.len, self->reference.str,
                    variant_from_index(self, node_idx), self->intervals.nodes[node_idx].distance,
                    supremal, distance);
                gva_string_destroy(allocator, supremal.sequence);

                hits[i].relation = relation_from_included(hits[i].included, self->intervals.nodes[node_idx].distance, distance, &hits[i].excluded);
                hits[i].query.start = start_ranges;
                hits[i].query.end = array_length(ranges);

                // disable aggregated hits except the last
                if (prev != GVA_NULL)
                {
                    hits[prev].next = i;
                } // if
                else
                {
                    entries[idx].head = i;
                } // else
                prev = i;

                entries[idx].included += hits[i].included;
                continue;
            } // if

            // check single query hit with multiple join hits
            while (hits[i].next != GVA_NULL && hits[i].query.start == hits[hits[i].next].query.start)
            {
                if (distance == 0)
                {
                    ARRAY_APPEND(allocator, ranges, hits[i].join.start);
                    gva_uint const node_idx = self->join[hits[i].join.start].link ^ entries[idx].gva_key;
                    supremal = gva_variant_dup(allocator, variant_from_index(self, node_idx));
                    distance = self->intervals.nodes[node_idx].distance;
                } // if
                i = hits[i].next;
                ARRAY_APPEND(allocator, ranges, hits[i].join.start);
                gva_uint const node_idx = self->join[hits[i].join.start].link ^ entries[idx].gva_key;
                GVA_Variant const variant = variant_from_index(self, node_idx);
                supremal.sequence = gva_string_concat(allocator, supremal.sequence,
                    (GVA_String) {variant.start - supremal.end, self->reference.str + supremal.end});
                supremal.sequence = gva_string_concat(allocator, supremal.sequence, variant.sequence);
                supremal.end = variant.end;
                distance += self->intervals.nodes[node_idx].distance;
                //fprintf(stderr, "  %2u  \" %2u (join)\n", hits[i].join.start, hits[i].included);
            } // while
            if (start_ranges < array_length(ranges))
            {
                //fprintf(stderr, "  SOLVE multiple join: " GVA_VARIANT_FMT_SPDI " (%zu)\n", GVA_VARIANT_PRINT_SPDI("NC_000001.11", supremal), distance);
                hits[i].included = variants_with_distance(allocator, self->reference.len, self->reference.str,
                    supremal, distance,
                    gva_lcs_graph_local_supremal(graph, hits[i].query.start, hits[i].query.end), graph.dom_nodes[hits[i].query.end].distance - graph.dom_nodes[hits[i].query.start].distance);
                gva_string_destroy(allocator, supremal.sequence);

                hits[i].relation = relation_from_included(hits[i].included, distance, graph.dom_nodes[hits[i].query.end].distance - graph.dom_nodes[hits[i].query.start].distance, &hits[i].excluded);
                hits[i].join.start = start_ranges;
                hits[i].join.end = array_length(ranges);

                // disable aggregated hits except the last
                if (prev != GVA_NULL)
                {
                    hits[prev].next = i;
                } // if
                else
                {
                    entries[idx].head = i;
                } // else
                prev = i;

                entries[idx].included += hits[i].included;
                continue;
            } // if

            prev = i;
            distance = graph.dom_nodes[hits[i].query.end].distance - graph.dom_nodes[hits[i].query.start].distance;

            gva_uint const node_idx = self->join[hits[i].join.start].link ^ entries[idx].gva_key;
            // non-distance based
            if (ABS((intmax_t) self->intervals.nodes[node_idx].distance - (intmax_t) distance) != hits[i].included)
            {
                //fprintf(stderr, "  SOLVE non-distance based\n");
                size_t excluded = 0;
                hits[i].included = variants_included(allocator, self->reference.len, self->reference.str,
                    variant_from_index(self, node_idx),
                    gva_lcs_graph_local_supremal(graph, hits[i].query.start, hits[i].query.end),
                    &excluded);
                if (hits[i].included > 0)
                {
                    hits[i].included = MAX(1, (double) (2 * hits[i].included - 1) / excluded * MIN(self->intervals.nodes[node_idx].distance, distance));
                } // if
                hits[i].relation = relation_from_included(hits[i].included, self->intervals.nodes[node_idx].distance, distance, &hits[i].excluded);

                entries[idx].included += hits[i].included;
                continue;
            } // if

            // equivalent or contains
            if (self->intervals.nodes[node_idx].distance - distance == hits[i].included)
            {
                hits[i].included = distance;
            } // if
            // is_contained
            else if (distance - self->intervals.nodes[node_idx].distance == hits[i].included)
            {
                hits[i].included = self->intervals.nodes[node_idx].distance;
            } // if
            hits[i].relation = relation_from_included(hits[i].included, self->intervals.nodes[node_idx].distance, distance, &hits[i].excluded);

            entries[idx].included += hits[i].included;
        } // for

        if (entries[idx].included == 0)
        {
            entries[idx].gva_key = NOT_FOUND;
            continue;
        } // if

        // disjoint alleles
        gva_uint* join = NULL;
        gva_uint prev_join = self->alleles[entries[idx].gva_key].join.start;
        size_t distance_join = 0;
        gva_uint* query = NULL;
        gva_uint prev_query = 0;
        size_t distance_query = 0;

        gva_uint excluded = 0;
        GVA_Relation const relation = relation_from_included(entries[idx].included, self->alleles[entries[idx].gva_key].distance, gva_lcs_graph_distance(graph), &excluded);
        fprintf(stderr, GVA_STRING_FMT " %2u %2u %s\n", GVA_STRING_PRINT(gva_index_id(self, entries[idx].gva_key)), entries[idx].included, excluded, GVA_RELATION_LABELS[relation]);
        for (gva_uint i = entries[idx].head; i != GVA_NULL; i = hits[i].next)
        {
            if (hits[i].join.end - hits[i].join.start > 1)
            {
                fprintf(stderr, "  [");
                for (gva_uint j = hits[i].join.start; j < hits[i].join.end; ++j)
                {
                    for (gva_uint k = prev_join; k < ranges[j]; ++k)
                    {
                        ARRAY_APPEND(allocator, join, k);
                        distance_join += self->intervals.nodes[self->join[k].link ^ entries[idx].gva_key].distance;
                    } // for
                    prev_join = ranges[j] + 1;
                    fprintf(stderr, "%2u, ", ranges[j]);
                } // for
                fprintf(stderr, "] %2u %2u %2u %s\n", hits[i].query.start, hits[i].included, hits[i].excluded, GVA_RELATION_LABELS[hits[i].relation]);
                for (gva_uint k = prev_query; k < hits[i].query.start; ++k)
                {
                    ARRAY_APPEND(allocator, query, k);
                } // for
                distance_query += graph.dom_nodes[hits[i].query.start].distance - graph.dom_nodes[prev_query].distance;
                prev_query = hits[i].query.start + 1;
            } // if
            else if (hits[i].query.end - hits[i].query.start > 1)
            {
                fprintf(stderr, "  %2u [", hits[i].join.start);
                for (gva_uint j = hits[i].query.start; j < hits[i].query.end; ++j)
                {
                    for (gva_uint k = prev_query; k < ranges[j]; ++k)
                    {
                        ARRAY_APPEND(allocator, query, k);
                    } // for
                    distance_query += graph.dom_nodes[ranges[j]].distance - graph.dom_nodes[prev_query].distance;
                    prev_query = ranges[j] + 1;
                    fprintf(stderr, "%2u, ", ranges[j]);
                } // for
                fprintf(stderr, "] %2u %2u %s\n", hits[i].included, hits[i].excluded, GVA_RELATION_LABELS[hits[i].relation]);
                for (gva_uint k = prev_join; k < hits[i].join.start; ++k)
                {
                    ARRAY_APPEND(allocator, join, k);
                    distance_join += self->intervals.nodes[self->join[k].link ^ entries[idx].gva_key].distance;
                } // for
                prev_join = hits[i].join.start + 1;
            } // if
            else if (hits[i].included > 0)
            {
                fprintf(stderr, "  %2u %2u %2u %2u %s\n", hits[i].join.start, hits[i].query.start, hits[i].included, hits[i].excluded, GVA_RELATION_LABELS[hits[i].relation]);
                for (gva_uint k = prev_join; k < hits[i].join.start; ++k)
                {
                    ARRAY_APPEND(allocator, join, k);
                    distance_join += self->intervals.nodes[self->join[k].link ^ entries[idx].gva_key].distance;
                } // for
                prev_join = hits[i].join.start + 1;
                for (gva_uint k = prev_query; k < hits[i].query.start; ++k)
                {
                    ARRAY_APPEND(allocator, query, k);
                } // for
                distance_query += graph.dom_nodes[hits[i].query.start].distance - graph.dom_nodes[prev_query].distance;
                prev_query = hits[i].query.start + 1;
            } // if
        } // for
        for (gva_uint k = prev_join; k < self->alleles[entries[idx].gva_key].join.end; ++k)
        {
            ARRAY_APPEND(allocator, join, k);
            distance_join += self->intervals.nodes[self->join[k].link ^ entries[idx].gva_key].distance;
        } // for
        for (gva_uint k = prev_query; k < array_length(graph.dom_nodes) - 1; ++k)
        {
            ARRAY_APPEND(allocator, query, k);
        } // for
        distance_query += graph.dom_nodes[array_length(graph.dom_nodes) - 1].distance - graph.dom_nodes[prev_query].distance;

        if (distance_join > 0)
        {
            fprintf(stderr, "  [");
            for (size_t i = 0; i < array_length(join); ++i)
            {
                fprintf(stderr, "%2u, ", join[i]);
            } // for
            fprintf(stderr, "]  *  0 %2zu disjoint\n", distance_join);
        } // if
        if (distance_query > 0)
        {
            fprintf(stderr, "   * [");
            for (size_t i = 0; i < array_length(query); ++i)
            {
                fprintf(stderr, "%2u, ", query[i]);
            } // for
            fprintf(stderr, "]  0 %2zu disjoint\n", distance_query);
        } // if

        join = ARRAY_DESTROY(allocator, join);
        query = ARRAY_DESTROY(allocator, query);
    } // for

    ranges = ARRAY_DESTROY(allocator, ranges);
    hits = ARRAY_DESTROY(allocator, hits);
    entries = HASH_TABLE_DESTROY(allocator, entries);

    return (GVA_Query_Result) {NULL};
} // gva_index_query


inline GVA_String
gva_index_id(GVA_Index const* const self, size_t const idx)
{
    return trie_string(self->ids, self->alleles[idx].id_idx);
} // gva_index_id


inline GVA_Variant
gva_index_variant(GVA_Index const* const self,
    size_t const allele_idx, size_t const variant_idx)
{
    return variant_from_index(self, self->join[variant_idx].link ^ allele_idx);
} // gva_index_variant


gva_uint
gva_index_allele_start(GVA_Index const* const self, size_t const allele_idx)
{
    return self->alleles[allele_idx].join.start;
} // gva_index_allele_start


gva_uint
gva_index_node_distance(GVA_Index const* const self, size_t const allele_idx, size_t const variant_idx)
{
    return self->intervals.nodes[self->join[variant_idx].link ^ allele_idx].distance;
} // gva_index_node_distance
