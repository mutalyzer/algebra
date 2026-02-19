#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // intmax_t
#include <string.h>     // memcpy

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/edit.h"        // gva_edit_distance
#include "../include/index.h"       // GVA_Index, gva_index_*, GVA_Result
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_lcs_graph_*
#include "../include/relations.h"   // GVA_Relation, GVA_CONTAINS,
                                    // GVA_DISJOINT, GVA_EQUIVALENT,
                                    // GVA_IS_CONTAINED
#include "../include/string.h"      // GVA_String
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "../include/variant.h"     // GVA_Variant, gva_variant_*
#include "array.h"              // ARRAY_*, array_*
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
    gva_uint id_idx;
    gva_uint distance;

    // wrt to join
    gva_uint start;
    gva_uint end;
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
    GVA_Variant const lhs, GVA_Variant const rhs)
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
    gva_lcs_graph_destroy(allocator, lhs_graph);

    gva_lcs_graph_uniq_atomics(rhs_graph, start, start_intersection, end_intersection, rhs_dels, rhs_as, rhs_cs, rhs_gs, rhs_ts);
    gva_lcs_graph_destroy(allocator, rhs_graph);

    size_t const included =
        bitset_intersection_cnt(lhs_dels, rhs_dels) +
        bitset_intersection_cnt(lhs_as, rhs_as) +
        bitset_intersection_cnt(lhs_cs, rhs_cs) +
        bitset_intersection_cnt(lhs_gs, rhs_gs) +
        bitset_intersection_cnt(lhs_ts, rhs_ts);

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

    return included > 0 ? 1 : 0;
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

    return variants_included(allocator, len_ref, reference, lhs, rhs);
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

    index->intervals = interval_tree_init();
    index->inserted = trie_init();
    index->ids = trie_init();
    index->join = NULL;
    index->alleles = NULL;

    return index;
} // gva_index_init


inline GVA_Index*
gva_index_destroy(GVA_Index* const self)
{
    if (self == NULL)
    {
        return NULL;
    } // if

    interval_tree_destroy(self->allocator, &self->intervals);
    trie_destroy(self->allocator, &self->inserted);
    trie_destroy(self->allocator, &self->ids);
    self->join = ARRAY_DESTROY(self->allocator, self->join);
    self->alleles = ARRAY_DESTROY(self->allocator, self->alleles);

    return self->allocator.allocate(self->allocator.context, self, sizeof(*self), 0);
} // gva_index_destroy


void
gva_index_insert(GVA_Index* restrict const self,
    size_t const len_id, char const id[static restrict len_id],
    GVA_Variant const variant, size_t const distance)
{
    size_t const id_idx = trie_insert(self->allocator, &self->ids, len_id, id);
    size_t allele_idx = array_length(self->alleles) - 1;

    // new allele
    if (array_length(self->alleles) == 0 || id_idx != self->alleles[allele_idx].id_idx)
    {
        allele_idx = ARRAY_APPEND(self->allocator, self->alleles, ((Allele) {id_idx, 0, array_length(self->join), 0}));
    } // if

    // add variant
    gva_uint const inserted_idx = trie_insert(self->allocator, &self->inserted, variant.sequence.len, variant.sequence.str);
    gva_uint const tmp_idx = ARRAY_APPEND(self->allocator, self->intervals.nodes,
                                          ((Interval_Tree_Node) {{GVA_NULL, GVA_NULL},
                                          variant.start, variant.end, variant.end, 0, inserted_idx, GVA_NULL, distance}));
    gva_uint const node_idx = interval_tree_insert(&self->intervals, tmp_idx);
    // undo append: interval already in the tree
    if (node_idx != tmp_idx)
    {
        array_header(self->intervals.nodes)->length -= 1;
    } // if
    self->intervals.nodes[node_idx].alleles = ARRAY_APPEND(self->allocator, self->join,
                                  ((Join) {node_idx ^ allele_idx, self->intervals.nodes[node_idx].alleles}));

    // update allele
    self->alleles[allele_idx].distance += distance;
    self->alleles[allele_idx].end = self->intervals.nodes[node_idx].alleles;
} // gva_index_insert


GVA_Result*
gva_index_query(GVA_Allocator const allocator,
    GVA_Index* const self, GVA_LCS_Graph const graph)
{
    struct Allele_Entry
    {
        HASH_TABLE_KEY;     // index into self->alleles
        gva_uint included;
        gva_uint head;      // singly linked list in hits
        gva_uint tail;
    }* entries = hash_table_init(allocator, 1024, sizeof(*entries));

    struct Hit
    {
        gva_uint idx;       // index into self->join
        gva_uint included;  // double purpose: distance or included
        gva_uint start;     // wrt local supremal parts in the query
        gva_uint end;
        gva_uint next;
    }* hits = NULL;

    //
    // Phase 1: for every local supremal part:
    //          - find candidates based on interval query
    //          - calculate the distance with the candidates,
    //            and quickly eject disjoint candidates based on distance
    //
    for (size_t i = 0; i < array_length(graph.dom_nodes) - 1; ++i)
    {
        GVA_Variant const variant = gva_lcs_graph_ls_slice(graph, i, i + 1);
        // fprintf(stderr, "    %2zu " GVA_VARIANT_FMT_SPDI " (%u)\n", i, GVA_VARIANT_PRINT_SPDI("NC_000001.11", variant), graph.dom_nodes[i + 1].distance);

        gva_uint* intervals = interval_tree_intersection(allocator, self->intervals, variant.start, variant.end);
        for (size_t j = 0; j < array_length(intervals); ++j)
        {
            gva_uint const node_idx = intervals[j];
            size_t const distance = variants_distance(allocator, self->reference.len, self->reference.str,
                                                      variant_from_index(self, node_idx), variant);
            if (distance >= self->intervals.nodes[node_idx].distance + graph.dom_nodes[i + 1].distance)
            {
                // Disjoint based on distance
                continue;
            } // if

            // fprintf(stderr, "        vs " GVA_VARIANT_FMT_SPDI " (%u)\n", GVA_VARIANT_PRINT_SPDI("NC_000001.11", variant_from_index(self, node_idx)), self->intervals.nodes[node_idx].distance);

            for (gva_uint k = self->intervals.nodes[node_idx].alleles; k != GVA_NULL; k = self->join[k].next)
            {
                gva_uint const allele_idx = self->join[k].link ^ node_idx;

                size_t const idx = HASH_TABLE_INDEX(entries, allele_idx);
                if (entries[idx].gva_key != allele_idx)
                {
                    // The first time a part is associated with this allele_idx
                    // we point it to the freshly created result part entry
                    gva_uint const tail = ARRAY_APPEND(allocator, hits,
                        ((struct Hit)
                        {
                            .idx = k,
                            .included = distance,
                            .start = i,
                            .end = i + 1,
                            .next = GVA_NULL
                        }));
                    HASH_TABLE_SET(allocator, entries, allele_idx,
                        ((struct Allele_Entry)
                        {
                            .gva_key = allele_idx,
                            .head = tail,
                            .tail = tail
                        }));
                    continue;
                } // if

                // is there a new hit for this allele?
                if (hits[entries[idx].tail].idx != k)
                {
                    gva_uint const tail = ARRAY_APPEND(allocator, hits,
                        ((struct Hit)
                        {
                            .idx = k,
                            .included = distance,
                            .start = i,
                            .end = i + 1,
                            .next = GVA_NULL
                        }));
                    // add to end of list
                    hits[entries[idx].tail].next = tail;
                    entries[idx].tail = tail;
                } // if
                // otherwise the same hit for this query part
                else
                {
                    hits[entries[idx].tail].end = i + 1;
                } // else
            } // for
        } // for
        intervals = ARRAY_DESTROY(allocator, intervals);
    } // for

    //
    // Phase 2:
    //
    for (size_t idx = 0; idx < array_header(entries)->capacity; ++idx)
    {
        if (entries[idx].gva_key == NOT_FOUND)
        {
            continue;
        } // if

        // fprintf(stderr, GVA_STRING_FMT ":\n", GVA_STRING_PRINT(trie_string(self->ids, self->alleles[entries[idx].gva_key].id_idx)));

        gva_uint prev = GVA_NULL;
        for (gva_uint i = entries[idx].head; i != GVA_NULL; i = hits[i].next)
        {
            // fprintf(stderr, "    %u %u  %u %u\n", hits[i].idx, hits[i].included, hits[i].start, hits[i].end);

            GVA_Variant lhs = {0, 0, {0, NULL}};
            size_t distance = 0;
            // look ahead for the same local supremal part in the query
            while (hits[i].next != GVA_NULL && hits[hits[i].next].start == hits[i].start)
            {
                if (gva_variant_length(lhs) == 0)
                {
                    gva_uint const node_idx = self->join[hits[i].idx].link ^ entries[idx].gva_key;
                    lhs = gva_variant_dup(allocator, variant_from_index(self, node_idx));
                    distance = self->intervals.nodes[node_idx].distance;
                } // if
                i = hits[i].next;
                gva_uint const node_idx = self->join[hits[i].idx].link ^ entries[idx].gva_key;
                GVA_Variant const variant = variant_from_index(self, node_idx);
                lhs.sequence = gva_string_concat(allocator, lhs.sequence,
                    (GVA_String) {variant.start - lhs.end, self->reference.str + lhs.end});
                lhs.sequence = gva_string_concat(allocator, lhs.sequence, variant.sequence);
                lhs.end = variant.end;
                distance += self->intervals.nodes[node_idx].distance;

                // fprintf(stderr, "    %u %u  %u %u\n", hits[i].idx, hits[i].included, hits[i].start, hits[i].end);
            } // while

            gva_uint const node_idx = self->join[hits[i].idx].link ^ entries[idx].gva_key;

            // single query part has multiple hits
            if (gva_variant_length(lhs) != 0)
            {
                // aggregate hits in last hit
                hits[i].included = variants_with_distance(allocator, self->reference.len, self->reference.str,
                    lhs, distance,
                    gva_lcs_graph_ls_slice(graph, hits[i].start, hits[i].end), graph.dom_nodes[hits[i].start + 1].distance);

                // disable aggregated hits except the last
                if (prev != GVA_NULL)
                {
                    hits[prev].next = i;
                } // if
                else
                {
                    entries[idx].head = i;
                } // else

                // fprintf(stderr, "        DB multi hit: " GVA_VARIANT_FMT_SPDI " (%zu) :: %u\n", GVA_VARIANT_PRINT_SPDI("NC_000001.11", lhs), distance, hits[i].included);
                gva_string_destroy(allocator, lhs.sequence);
            } // if
            // single hit with multiple query parts
            else if (hits[i].end - hits[i].start > 1)  // FIXME: similar body to if below?
            {
                // FIXME: use cumulative distances in dom_nodes
                size_t distance = 0;  // graph.dom_nodes[hits[i].end].distance - graph.dom_nodes[hits[i].start].distance
                for (gva_uint j = hits[i].start; j < hits[i].end; ++j)
                {
                    distance += graph.dom_nodes[j + 1].distance;
                } // for
                hits[i].included = variants_with_distance(allocator, self->reference.len, self->reference.str,
                    variant_from_index(self, node_idx), self->intervals.nodes[node_idx].distance,
                    gva_lcs_graph_ls_slice(graph, hits[i].start, hits[i].end), distance);

                // fprintf(stderr, "        Query multi hit: %u (%zu) :: %u\n", hits[i].end - hits[i].start, distance, hits[i].included);
            } // if
            // single hit for single query part, relation is not determined by the distances
            else if (ABS((intmax_t) self->intervals.nodes[node_idx].distance - graph.dom_nodes[hits[i].start + 1].distance) != hits[i].included)
            {
                hits[i].included = variants_included(allocator, self->reference.len, self->reference.str,
                    variant_from_index(self, node_idx),
                    gva_lcs_graph_ls_slice(graph, hits[i].start, hits[i].end));
                // fprintf(stderr, "        Fully calculate: %u\n", hits[i].included);
            } // if
            // single hit for single query part, relation is already determined
            else
            {
                // equivalent or contains
                if (self->intervals.nodes[node_idx].distance - graph.dom_nodes[hits[i].start + 1].distance == hits[i].included)
                {
                    hits[i].included = graph.dom_nodes[hits[i].start + 1].distance;
                } // if
                // is_contained
                else if (graph.dom_nodes[hits[i].start + 1].distance - self->intervals.nodes[node_idx].distance == hits[i].included)
                {
                    hits[i].included = self->intervals.nodes[node_idx].distance;
                } // if
            } // else
            entries[idx].included += hits[i].included;
            prev = i;
        } // for
        // TODO: exclude included == 0?
        if (entries[idx].included == 0)
        {
            entries[idx].gva_key = NOT_FOUND;
        } // if
    } // for

    //
    // Phase 3:
    //

    GVA_Result* results = NULL;

    for (size_t idx = 0; idx < array_header(entries)->capacity; ++idx)
    {
        if (entries[idx].gva_key == NOT_FOUND)
        {
            continue;
        } // if


        // fprintf(stderr, GVA_STRING_FMT ": %u\n",
        //         GVA_STRING_PRINT(trie_string(self->ids, self->alleles[entries[idx].gva_key].id_idx)),
        //         entries[idx].included);

        // for (gva_uint i = entries[idx].head; i != GVA_NULL; i = hits[i].next)
        // {
        //     GVA_Variant variant;
        //     gva_edges(graph.observed.str,
        //                 graph.dom_nodes[hits[i].start], graph.dom_nodes[hits[i].end],
        //                 hits[i].start == 0, hits[i].end == array_length(graph.dom_nodes) - 2,
        //                 &variant);
        //
        //     gva_uint const node_idx = self->join[hits[i].idx].link ^ entries[idx].gva_key;
        //     fprintf(stderr, "    " GVA_VARIANT_FMT_SPDI " " GVA_VARIANT_FMT_SPDI " %u\n", GVA_VARIANT_PRINT_SPDI("NC_000006.12", variant_from_index(self, node_idx)), GVA_VARIANT_PRINT_SPDI("NC_000006.12", variant), hits[i].included);
        // } // for

        GVA_Relation relation;
        size_t const excluded = self->alleles[entries[idx].gva_key].distance + graph.distance - 2 * entries[idx].included;
        if (excluded == 0)
        {
            relation = GVA_EQUIVALENT;
        } // if
        else if (entries[idx].included == graph.distance)
        {
            relation = GVA_CONTAINS;
        } // else if
        else if (entries[idx].included == self->alleles[entries[idx].gva_key].distance)
        {
            relation = GVA_IS_CONTAINED;
        } // else if
        else
        {
            relation = GVA_OVERLAP;
        } // else

        ARRAY_APPEND(self->allocator, results,
            ((GVA_Result) {trie_string(self->ids, self->alleles[entries[idx].gva_key].id_idx), relation, entries[idx].included, excluded})
        );
    } // for

    hits = ARRAY_DESTROY(allocator, hits);
    entries = HASH_TABLE_DESTROY(allocator, entries);

    return results;
} // gva_index_query
