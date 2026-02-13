#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // intmax_t
#include <string.h>     // memcpy


#include "../include/allocator.h"   // GVA_Allocator
#include "../include/edit.h"        // gva_edit_distance
#include "../include/index.h"       // GVA_Index, gva_index_*
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_edges
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
    memcpy((char*) observed_lhs.str + lhs.start - start, lhs.sequence.str, lhs.sequence.len);
    memcpy((char*) observed_lhs.str + lhs.start - start + lhs.sequence.len, reference + lhs.end, end - lhs.end);

    memcpy((char*) observed_rhs.str, reference + start, rhs.start - start);
    memcpy((char*) observed_rhs.str + rhs.start - start, rhs.sequence.str, rhs.sequence.len);
    memcpy((char*) observed_rhs.str + rhs.start - start + rhs.sequence.len, reference + rhs.end, end - rhs.end);

    size_t const distance = gva_edit_distance(allocator, observed_lhs.len, observed_lhs.str, observed_rhs.len, observed_rhs.str);

    gva_string_destroy(allocator, observed_rhs);
    gva_string_destroy(allocator, observed_lhs);

    return distance;
} // variants_distance


// FIXME: taken from compare_supremals
// FIXME: overlap with variants_distance
static size_t
variants_common(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref],
    GVA_Variant const lhs, GVA_Variant const rhs)
{
    size_t const start = MIN(lhs.start, rhs.start);
    size_t const end = MAX(lhs.end, rhs.end);

    GVA_String observed_lhs = gva_string_init(allocator, (lhs.start - start) + lhs.sequence.len + (end - lhs.end));
    GVA_String observed_rhs = gva_string_init(allocator, (rhs.start - start) + rhs.sequence.len + (end - rhs.end));

    memcpy((char*) observed_lhs.str, reference + start, lhs.start - start);
    memcpy((char*) observed_lhs.str + lhs.start - start, lhs.sequence.str, lhs.sequence.len);
    memcpy((char*) observed_lhs.str + lhs.start - start + lhs.sequence.len, reference + lhs.end, end - lhs.end);

    memcpy((char*) observed_rhs.str, reference + start, rhs.start - start);
    memcpy((char*) observed_rhs.str + rhs.start - start, rhs.sequence.str, rhs.sequence.len);
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

    // could be done on intersection instead of union
    bitset_fill(rhs_graph, start, start_intersection, end_intersection, rhs_dels, rhs_as, rhs_cs, rhs_gs, rhs_ts);
    gva_lcs_graph_destroy(allocator, rhs_graph);

    bitset_fill(lhs_graph, start, start_intersection, end_intersection, lhs_dels, lhs_as, lhs_cs, lhs_gs, lhs_ts);
    gva_lcs_graph_destroy(allocator, lhs_graph);

    size_t const common =
        bitset_intersection_cnt(lhs_dels, rhs_dels) +
        bitset_intersection_cnt(lhs_as, rhs_as) +
        bitset_intersection_cnt(lhs_cs, rhs_cs) +
        bitset_intersection_cnt(lhs_gs, rhs_gs) +
        bitset_intersection_cnt(lhs_ts, rhs_ts);

    // Also calculate the union
    // Note: bitset_fill needs to use "plain" start and end.
    //
    // size_t const union1 =
    //         bitset_intersection_cnt(lhs_dels, lhs_dels) +
    //         bitset_intersection_cnt(lhs_as, lhs_as) +
    //         bitset_intersection_cnt(lhs_cs, lhs_cs) +
    //         bitset_intersection_cnt(lhs_gs, lhs_gs) +
    //         bitset_intersection_cnt(lhs_ts, lhs_ts) +
    //         bitset_intersection_cnt(rhs_dels, rhs_dels) +
    //         bitset_intersection_cnt(rhs_as, rhs_as) +
    //         bitset_intersection_cnt(rhs_cs, rhs_cs) +
    //         bitset_intersection_cnt(rhs_gs, rhs_gs) +
    //         bitset_intersection_cnt(rhs_ts, rhs_ts) - common;

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

    return common > 0 ? 1 : 0;
} // variants_common


inline GVA_Index*
gva_index_init(GVA_Allocator const allocator,
    size_t const len, char const reference[static len])
{
    GVA_Index* index = allocator.allocate(allocator.context, NULL, 0, sizeof(*index));
    if (index == NULL)
    {
        return NULL;
    } // if

    index->allocator = allocator;
    index->reference = (GVA_String) {len, reference};

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
    size_t const len, char const allele_id[static restrict len],
    GVA_Variant const variant, size_t const distance)
{
    size_t const id_idx = trie_insert(self->allocator, &self->ids, len, allele_id);
    size_t allele_idx = array_length(self->alleles) - 1;

    // new allele
    if (array_length(self->alleles) == 0 || id_idx != self->alleles[allele_idx].id_idx)
    {
        allele_idx = ARRAY_APPEND(self->allocator, self->alleles, ((Allele) {id_idx, 0, array_length(self->join), 0})) - 1;
    } // if

    // add variant
    gva_uint const inserted_idx = trie_insert(self->allocator, &self->inserted, variant.sequence.len, variant.sequence.str);
    gva_uint const tmp_idx = ARRAY_APPEND(self->allocator, self->intervals.nodes,
                                          ((Interval_Tree_Node) {{GVA_NULL, GVA_NULL},
                                          variant.start, variant.end, variant.end, 0, inserted_idx, GVA_NULL, distance})) - 1;
    gva_uint const node_idx = interval_tree_insert(&self->intervals, tmp_idx);
    // undo append: interval already in the tree
    if (node_idx != tmp_idx)
    {
        array_header(self->intervals.nodes)->length -= 1;
    } // if
    self->intervals.nodes[node_idx].alleles = ARRAY_APPEND(self->allocator, self->join, ((Join) {node_idx ^ allele_idx, self->intervals.nodes[node_idx].alleles})) - 1;

    // update allele
    self->alleles[allele_idx].distance += distance;
    self->alleles[allele_idx].end = self->intervals.nodes[node_idx].alleles;
} // gva_index_insert


void
gva_index_query(GVA_Allocator const allocator,
    GVA_Index* const self, GVA_LCS_Graph const graph)
{
    struct Result_Allele
    {
        HASH_TABLE_KEY;
        gva_uint included;
        gva_uint head;
        gva_uint tail;
    }* alleles = hash_table_init(allocator, 1024, sizeof(*alleles));

    struct Result_Part
    {
        gva_uint part;
        gva_uint distance;
        gva_uint start;
        gva_uint end;
        gva_uint next;
    }* parts = NULL;

    for (size_t i = 0; i < array_length(graph.dom_nodes) - 1; ++i)
    {
        GVA_Variant variant;
        gva_edges(graph.observed.str,
                  graph.dom_nodes[i], graph.dom_nodes[i + 1],
                  i == 0, i == array_length(graph.dom_nodes) - 2,
                  &variant);

        fprintf(stderr, "    %2zu " GVA_VARIANT_FMT_SPDI " (%u)\n", i, GVA_VARIANT_PRINT_SPDI("NC_000001.11", variant), graph.dom_nodes[i + 1].distance);

        gva_uint* intervals = interval_tree_intersection(allocator, self->intervals, variant.start, variant.end);
        for (size_t j = 0; j < array_length(intervals); ++j)
        {
            size_t const distance = variants_distance(allocator, self->reference.len, self->reference.str, variant_from_index(self, intervals[j]), variant);
            if (distance >= self->intervals.nodes[intervals[j]].distance + graph.dom_nodes[i + 1].distance)
            {
                continue;  // disjoint
            } // if

            fprintf(stderr, "        vs " GVA_VARIANT_FMT_SPDI " (%u)\n", GVA_VARIANT_PRINT_SPDI("NC_000001.11", variant_from_index(self, intervals[j])), self->intervals.nodes[intervals[j]].distance);

            for (gva_uint k = self->intervals.nodes[intervals[j]].alleles; k != GVA_NULL; k = self->join[k].next)
            {
                gva_uint const allele_idx = self->join[k].link ^ intervals[j];

                size_t const idx = HASH_TABLE_INDEX(alleles, allele_idx);
                if (alleles[idx].gva_key != allele_idx)
                {
                    gva_uint const tail = ARRAY_APPEND(allocator, parts, ((struct Result_Part) {k, distance, i, i + 1, GVA_NULL})) - 1;
                    HASH_TABLE_SET(allocator, alleles, allele_idx, ((struct Result_Allele) {allele_idx, 0, tail, tail}));
                    continue;
                } // if

                if (parts[alleles[idx].tail].part != k || parts[alleles[idx].tail].start == i)
                {
                    gva_uint const tail = ARRAY_APPEND(allocator, parts, ((struct Result_Part) {k, distance, i, i + 1, GVA_NULL})) - 1;
                    parts[alleles[idx].tail].next = tail;
                    alleles[idx].tail = tail;
                } // if
                else
                {
                    parts[alleles[idx].tail].end = i + 1;
                } // else
            } // for
        } // for
        intervals = ARRAY_DESTROY(allocator, intervals);
    } // for

    for (size_t idx = 0; idx < array_header(alleles)->capacity; ++idx)
    {
        if (alleles[idx].gva_key == NOT_FOUND)
        {
            continue;
        } // if

        fprintf(stderr, GVA_STRING_FMT ":\n", GVA_STRING_PRINT(trie_string(self->ids, self->alleles[alleles[idx].gva_key].id_idx)));

        for (gva_uint i = alleles[idx].head; i != GVA_NULL; i = parts[i].next)
        {
            fprintf(stderr, "    %u %u  %u %u\n", parts[i].part, parts[i].distance, parts[i].start, parts[i].end);

            GVA_Variant* variants = NULL;
            size_t distance = 0;
            while (parts[i].next != GVA_NULL && parts[parts[i].next].start == parts[i].start)
            {
                if (variants == NULL)
                {
                    gva_uint const node_idx = self->join[parts[i].part].link ^ alleles[idx].gva_key;
                    ARRAY_APPEND(allocator, variants, variant_from_index(self, node_idx));
                    distance += self->intervals.nodes[node_idx].distance;
                } // if
                i = parts[i].next;
                gva_uint const node_idx = self->join[parts[i].part].link ^ alleles[idx].gva_key;
                ARRAY_APPEND(allocator, variants, variant_from_index(self, node_idx));
                distance += self->intervals.nodes[node_idx].distance;
                fprintf(stderr, "    %u %u  %u %u\n", parts[i].part, parts[i].distance, parts[i].start, parts[i].end);
            } // while
            if (variants != NULL)
            {
                // TODO: multiple DB hits with one query hit
                fprintf(stderr, "        DB multi hit: %zu (%zu)\n", array_length(variants), distance);
                variants = ARRAY_DESTROY(allocator, variants);
            } // if
            else if (parts[i].end - parts[i].start > 1)
            {
                size_t distance = 0;
                for (gva_uint j = parts[i].start; j < parts[i].end; ++j)
                {
                    distance += graph.dom_nodes[j + 1].distance;
                } // for
                // TODO: multiple query hits with one DB hit
                fprintf(stderr, "        Query multi hit: %u (%zu)\n", parts[i].end - parts[i].start, distance);
            } // if
            else if (ABS((intmax_t) self->intervals.nodes[self->join[parts[i].part].link ^ alleles[idx].gva_key].distance - graph.dom_nodes[parts[i].start + 1].distance) != parts[i].distance)
            {
                GVA_Variant variant;
                gva_edges(graph.observed.str,
                          graph.dom_nodes[parts[i].start], graph.dom_nodes[parts[i].end],
                          parts[i].start == 0, parts[i].end == array_length(graph.dom_nodes) - 2,
                          &variant);
                gva_uint const included = variants_common(allocator, self->reference.len, self->reference.str, variant_from_index(self, self->join[parts[i].part].link ^ alleles[idx].gva_key), variant);
                fprintf(stderr, "        Fully calculate: %u\n", included);
            } // if
            // TODO: aggregate independent results
        } // for
    } // for

    parts = ARRAY_DESTROY(allocator, parts);
    alleles = HASH_TABLE_DESTROY(allocator, alleles);
} // gva_index_query
