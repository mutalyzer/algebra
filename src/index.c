#include <stddef.h>     // NULL, size_t
#include <string.h>     // memcpy


#include "../include/allocator.h"   // GVA_Allocator
#include "../include/compare.h"     // gva_compare_with_distance
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
#include "common.h"             // MAX, MIN
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

    char* lhs_obs = allocator.allocate(allocator.context, NULL, 0, len_lhs);
    char* rhs_obs = allocator.allocate(allocator.context, NULL, 0, len_rhs);
    if (lhs_obs == NULL || rhs_obs == NULL)
    {
        rhs_obs = allocator.allocate(allocator.context, rhs_obs, len_rhs, 0);
        lhs_obs = allocator.allocate(allocator.context, lhs_obs, len_lhs, 0);
        return -1;  // FIXME: OOM
    } // if

    memcpy(lhs_obs, reference + start, lhs.start - start);
    memcpy(lhs_obs + lhs.start - start, lhs.sequence.str, lhs.sequence.len);
    memcpy(lhs_obs + lhs.start - start + lhs.sequence.len, reference + lhs.end, end - lhs.end);

    memcpy(rhs_obs, reference + start, rhs.start - start);
    memcpy(rhs_obs + rhs.start - start, rhs.sequence.str, rhs.sequence.len);
    memcpy(rhs_obs + rhs.start - start + rhs.sequence.len, reference + rhs.end, end - rhs.end);

    size_t const distance = gva_edit_distance(allocator, len_lhs, lhs_obs, len_rhs, rhs_obs);

    rhs_obs = allocator.allocate(allocator.context, rhs_obs, len_rhs, 0);
    lhs_obs = allocator.allocate(allocator.context, lhs_obs, len_lhs, 0);

    return distance;
} // variants_distance


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

        fprintf(stderr, "    " GVA_VARIANT_FMT_SPDI "\n", GVA_VARIANT_PRINT_SPDI("NC_000001.11", variant));

        gva_uint* intervals = interval_tree_intersection(allocator, self->intervals, variant.start, variant.end);
        for (size_t j = 0; j < array_length(intervals); ++j)
        {
            size_t const distance = variants_distance(allocator, self->reference.len, self->reference.str, variant_from_index(self, intervals[j]), variant);
            if (distance >= self->intervals.nodes[intervals[j]].distance + graph.dom_nodes[i + 1].distance)
            {
                continue;  // disjoint
            } // if

            fprintf(stderr, "        vs " GVA_VARIANT_FMT_SPDI "\n", GVA_VARIANT_PRINT_SPDI("NC_000001.11", variant_from_index(self, intervals[j])));

            for (gva_uint k = self->intervals.nodes[intervals[j]].alleles; k != GVA_NULL; k = self->join[k].next)
            {
                gva_uint const allele_idx = self->join[k].link ^ intervals[j];

                size_t const idx = HASH_TABLE_INDEX(alleles, allele_idx);
                if (alleles[idx].gva_key != allele_idx)
                {
                    gva_uint const tail = ARRAY_APPEND(allocator, parts, ((struct Result_Part) {k, distance, i, i + 1, GVA_NULL})) - 1;
                    HASH_TABLE_SET(allocator, alleles, allele_idx, ((struct Result_Allele) {allele_idx, tail, tail}));
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
        } // for
    } // for

    parts = ARRAY_DESTROY(allocator, parts);
    alleles = HASH_TABLE_DESTROY(allocator, alleles);
} // gva_index_query
