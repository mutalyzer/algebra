#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // intmax_t, uint8_t

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/compare.h"     // gva_compare_*
#include "../include/edit.h"        // gva_edit_distance
#include "../include/index.h"       // GVA_Index, GVA_Query_Result, gva_index_*
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_lcs_graph_*
#include "../include/relations.h"   // GVA_Relation, GVA_CONTAINS,
                                    // GVA_DISJOINT, GVA_EQUIVALENT,
                                    // GVA_IS_CONTAINED
#include "../include/string.h"      // GVA_String
#include "../include/types.h"       // GVA_NULL, GVA_Interval, gva_uint
#include "../include/variant.h"     // GVA_Variant, gva_patch, gva_variant_dup
#include "array.h"              // ARRAY_*, array_length
#include "common.h"             // ABS, MAX, MIN
#include "dfa.h"                // dfa_*
#include "hash_table.h"         // HASH_TABLE_*, hash_table_*
#include "interval_tree.h"      // Interval_Tree, interval_tree_*
#include "trie.h"               // Trie, trie_*


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
    Trie          dfas;
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


static inline uint8_t*
dfa_from_index(GVA_Index const self[static 1], size_t const idx)
{
    return (uint8_t*) trie_string(self->dfas, self->intervals.nodes[idx].dfa).str;
} // dfa_from_index


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

    index->intervals = interval_tree_init(allocator);
    index->inserted = trie_init(allocator);
    index->dfas = trie_init(allocator);

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
    trie_destroy(&self->dfas);
    trie_destroy(&self->ids);
    self->alleles = ARRAY_DESTROY(self->allocator, self->alleles);
    self->join = ARRAY_DESTROY(self->allocator, self->join);

    return self->allocator.allocate(self->allocator.context, self, sizeof(*self), 0);
} // gva_index_destroy


void
gva_index_insert(GVA_Index* restrict const self,
    size_t const len, char const id[static restrict len],
    GVA_Variant const variant, size_t const distance)
{
    size_t const id_idx = trie_insert(&self->ids, len, id);
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
    if (node_idx == tmp_idx)
    {
        uint8_t* dfa = dfa_from_alignment(self->allocator, NULL,
            variant.end - variant.start, self->reference.str + variant.start,
            variant.sequence.len, variant.sequence.str);

        self->intervals.nodes[node_idx].dfa = trie_insert(&self->dfas, array_length(dfa), (char*) dfa);

        dfa = ARRAY_DESTROY(self->allocator, dfa);
    } // if
    else
    {
        // undo append: interval already in the tree
        array_header(self->intervals.nodes)->length -= 1;
    } // else
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
    if (entries == NULL)
    {
        return (GVA_Query_Result) {NULL};  // OOM
    } // if

    struct Hit
    {
        gva_uint join;
        gva_uint query;
        gva_uint distance;
        gva_uint next;
    }* hits = NULL;

    gva_uint* restrict starts = allocator.allocate(allocator.context, NULL, 0, sizeof(*starts) * (array_length(graph.dom_nodes) - 1));
    if (starts == NULL)
    {
        entries = HASH_TABLE_DESTROY(allocator, entries);
        return (GVA_Query_Result) {NULL};  // OOM
    } // if

    uint8_t* restrict dfas = NULL;

    // Phase 1: for every local supremal part:
    //          - find candidates based on interval query
    //          - calculate the distance with the candidates,
    //            and discard disjoint candidates early based on distance
    for (size_t i = 0; i < array_length(graph.dom_nodes) - 1; ++i)
    {
        starts[i] = array_length(dfas);
        dfas = dfa_from_lcs_graph(allocator, dfas, graph, i, i + 1);

        GVA_Variant const variant = gva_lcs_graph_local_supremal(graph, i, i + 1);

        gva_uint* intervals = interval_tree_intersection(allocator, self->intervals, variant.start, variant.end);
        for (size_t j = 0; j < array_length(intervals); ++j)
        {
            gva_uint const node_idx = intervals[j];

            size_t const distance = gva_compare_distance(allocator, self->reference.len, self->reference.str,
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
                        .join = k,
                        .query = i,
                        .distance = distance,
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

    GVA_Query_Result result = {NULL};

    // Phase 2
    for (size_t idx = 0; idx < array_header(entries)->capacity; ++idx)
    {
        if (entries[idx].gva_key == HASH_TABLE_NOT_FOUND)
        {
            continue;
        } // if

        size_t const start_hits = array_length(result.hits);
        for (gva_uint i = entries[idx].head; i != GVA_NULL; i = hits[i].next)
        {
            // look aheads for multi hits: build new local supremal
            size_t const start_parts = array_length(result.parts);
            GVA_Variant supremal = {0};
            size_t distance = 0;

            // check single join hit with multiple query hits
            gva_uint start = hits[i].query;
            gva_uint end = hits[i].query + 1;
            while (hits[i].next != GVA_NULL && hits[i].join == hits[hits[i].next].join)
            {
                if (start_parts == array_length(result.parts))
                {
                    ARRAY_APPEND(allocator, result.parts, hits[i].query);
                } // if
                i = hits[i].next;
                ARRAY_APPEND(allocator, result.parts, hits[i].query);
                if (hits[i].query > end)
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
                    start = hits[i].query;
                } // if
                end = hits[i].query + 1;
            } // while
            if (start_parts < array_length(result.parts))
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
                gva_uint const node_idx = self->join[hits[i].join].link ^ entries[idx].gva_key;

                uint8_t* restrict dfa = dfa_init(allocator, supremal.end - supremal.start + 1, supremal.sequence.len + 1);
                size_t offset = 0;
                GVA_Variant prev = {0};
                for (size_t j = start_parts; j < array_length(result.parts); ++j)
                {
                    GVA_Variant const variant = gva_lcs_graph_local_supremal(graph, result.parts[j], result.parts[j] + 1);
                    if (j > start_parts)
                    {
                        offset += prev.sequence.len + variant.start - prev.end;
                    } // if
                    dfa = dfa_concat(supremal, dfa, variant, dfas + starts[result.parts[j]], offset);
                    prev = variant;
                } // for

                size_t const included = gva_compare_included(allocator, self->reference.len, self->reference.str,
                    variant_from_index(self, node_idx), self->intervals.nodes[node_idx].distance, dfa_from_index(self, node_idx),
                    supremal, distance, dfa);
                gva_string_destroy(allocator, supremal.sequence);

                size_t excluded = 0;
                GVA_Relation const relation = gva_compare_relation_from_included(included, self->intervals.nodes[node_idx].distance, distance, &excluded);
                ARRAY_APPEND(allocator, result.hits, ((GVA_Query_Hit)
                    {
                        .relation = relation,
                        .included = included,
                        .excluded = excluded,
                        .query_parts = {start_parts, array_length(result.parts)},
                        .index_parts = {hits[i].join, hits[i].join},
                    }));

                entries[idx].included += included;

                dfa = ARRAY_DESTROY(allocator, dfa);

                continue;
            } // if

            // check single query hit with multiple join hits
            while (hits[i].next != GVA_NULL && hits[i].query == hits[hits[i].next].query)
            {
                if (distance == 0)
                {
                    ARRAY_APPEND(allocator, result.parts, hits[i].join);
                    gva_uint const node_idx = self->join[hits[i].join].link ^ entries[idx].gva_key;
                    supremal = gva_variant_dup(allocator, variant_from_index(self, node_idx));
                    distance = self->intervals.nodes[node_idx].distance;
                } // if
                i = hits[i].next;
                ARRAY_APPEND(allocator, result.parts, hits[i].join);
                gva_uint const node_idx = self->join[hits[i].join].link ^ entries[idx].gva_key;
                GVA_Variant const variant = variant_from_index(self, node_idx);
                supremal.sequence = gva_string_concat(allocator, supremal.sequence,
                    (GVA_String) {variant.start - supremal.end, self->reference.str + supremal.end});
                supremal.sequence = gva_string_concat(allocator, supremal.sequence, variant.sequence);
                supremal.end = variant.end;

                distance += self->intervals.nodes[node_idx].distance;
            } // while
            if (start_parts < array_length(result.parts))
            {
                uint8_t* restrict dfa = dfa_init(allocator, supremal.end - supremal.start + 1, supremal.sequence.len + 1);
                size_t offset = 0;
                GVA_Variant prev = {0};
                for (size_t j = start_parts; j < array_length(result.parts); ++j)
                {
                    gva_uint const node_idx = self->join[result.parts[j]].link ^ entries[idx].gva_key;
                    GVA_Variant const variant = variant_from_index(self, node_idx);
                    if (j > start_parts)
                    {
                        offset += prev.sequence.len + variant.start - prev.end;
                    } // if
                    dfa = dfa_concat(supremal, dfa, variant, dfa_from_index(self, node_idx), offset);
                    prev = variant;
                } // for

                size_t const included = gva_compare_included(allocator, self->reference.len, self->reference.str,
                    supremal, distance, dfa,
                    gva_lcs_graph_local_supremal(graph, hits[i].query, hits[i].query + 1), graph.dom_nodes[hits[i].query + 1].distance - graph.dom_nodes[hits[i].query].distance, dfas + starts[hits[i].query]);
                gva_string_destroy(allocator, supremal.sequence);

                size_t excluded = 0;
                GVA_Relation const relation = gva_compare_relation_from_included(included, distance, graph.dom_nodes[hits[i].query + 1].distance - graph.dom_nodes[hits[i].query].distance, &excluded);

                ARRAY_APPEND(allocator, result.hits, ((GVA_Query_Hit)
                    {
                        .relation = relation,
                        .included = included,
                        .excluded = excluded,
                        .query_parts = {hits[i].query, hits[i].query},
                        .index_parts = {start_parts, array_length(result.parts)},
                    }));

                entries[idx].included += included;

                dfa = ARRAY_DESTROY(allocator, dfa);

                continue;
            } // if

            distance = graph.dom_nodes[hits[i].query + 1].distance - graph.dom_nodes[hits[i].query].distance;

            gva_uint const node_idx = self->join[hits[i].join].link ^ entries[idx].gva_key;
            // non-distance based
            if (ABS((intmax_t) self->intervals.nodes[node_idx].distance - (intmax_t) distance) != hits[i].distance)
            {
                //size_t const included = !dfa_disjoint(variant_from_index(self, node_idx), dfa_from_index(self, node_idx),
                //    gva_lcs_graph_local_supremal(graph, hits[i].query, hits[i].query + 1), dfas + starts[hits[i].query]);
                size_t const included = dfa_overlap(allocator, self->reference.len, self->reference.str,
                    variant_from_index(self, node_idx), dfa_from_index(self, node_idx),
                    gva_lcs_graph_local_supremal(graph, hits[i].query, hits[i].query + 1), dfas + starts[hits[i].query]);

                if (included > 0)
                {
                    ARRAY_APPEND(allocator, result.hits, ((GVA_Query_Hit)
                        {
                            .relation = GVA_OVERLAP,
                            .included = included,
                            .excluded = self->intervals.nodes[node_idx].distance + distance - 2 * included,
                            .query_parts = {hits[i].query, hits[i].query},
                            .index_parts = {hits[i].join, hits[i].join},
                        }));
                } // if

                entries[idx].included += included;
                continue;
            } // if

            size_t included = distance;  // equivalent or contains
            if (distance - self->intervals.nodes[node_idx].distance == hits[i].distance)
            {
                // is_contained
                included = self->intervals.nodes[node_idx].distance;
            } // if

            size_t excluded = 0;
            GVA_Relation const relation = gva_compare_relation_from_included(included, self->intervals.nodes[node_idx].distance, distance, &excluded);

            ARRAY_APPEND(allocator, result.hits, ((GVA_Query_Hit)
                {
                    .relation = relation,
                    .included = included,
                    .excluded = excluded,
                    .query_parts = {hits[i].query, hits[i].query},
                    .index_parts = {hits[i].join, hits[i].join},
                }));

            entries[idx].included += included;
        } // for

        if (entries[idx].included == 0)
        {
            continue;
        } // if

        size_t excluded = 0;
        GVA_Relation const relation = gva_compare_relation_from_included(entries[idx].included, self->alleles[entries[idx].gva_key].distance, gva_lcs_graph_distance(graph), &excluded);

        ARRAY_APPEND(allocator, result.alleles, ((GVA_Query_Allele)
            {
                .idx = entries[idx].gva_key,
                .relation = relation,
                .included = entries[idx].included,
                .excluded = excluded,
                .hits = {start_hits, array_length(result.hits)},
            }));
    } // for

    hits = ARRAY_DESTROY(allocator, hits);
    entries = HASH_TABLE_DESTROY(allocator, entries);
    dfas = ARRAY_DESTROY(allocator, dfas);
    starts = allocator.allocate(allocator.context, starts, sizeof(*starts) * (array_length(graph.dom_nodes) - 1), 0);

    return result;
} // gva_index_query


inline GVA_String
gva_index_allele_id(GVA_Index const* const self, size_t const idx)
{
    return trie_string(self->ids, self->alleles[idx].id_idx);
} // gva_index_allele_id


inline GVA_Interval
gva_index_allele_parts(GVA_Index const* const self, size_t const idx)
{
    return self->alleles[idx].join;
} // gva_index_allele_parts


inline GVA_Variant
gva_index_variant(GVA_Index const* const self,
    size_t const allele_idx, size_t const variant_idx)
{
    return variant_from_index(self, self->join[variant_idx].link ^ allele_idx);
} // gva_index_variant


inline size_t
gva_index_variant_distance(GVA_Index const* const self,
    size_t const allele_idx, size_t const variant_idx)
{
    return self->intervals.nodes[self->join[variant_idx].link ^ allele_idx].distance;
} // gva_index_variant_distance


void
gva_index_stats(GVA_Index* const self)
{
    fprintf(stderr, "reference: #%d, %zu\n", 1, self->reference.len);
    fprintf(stderr, "interval nodes: #%zu, %zu\n", array_length(self->intervals.nodes),
        array_length(self->intervals.nodes) * sizeof(Interval_Tree_Node));
    fprintf(stderr, "inserted nodes: #%zu, %zu\n", array_length(self->inserted.nodes),
        array_length(self->inserted.nodes) * sizeof(Trie_Node));
    fprintf(stderr, "inserted strings: #%zu, %zu\n", array_length(self->inserted.strings),
        array_length(self->inserted.strings) * sizeof(char));
    fprintf(stderr, "dfas nodes: #%zu, %zu\n", array_length(self->dfas.nodes),
        array_length(self->dfas.nodes) * sizeof(Trie_Node));
    fprintf(stderr, "dfas strings: #%zu, %zu\n", array_length(self->dfas.strings),
        array_length(self->dfas.strings) * sizeof(char));
    fprintf(stderr, "ids nodes: #%zu, %zu\n", array_length(self->ids.nodes),
        array_length(self->ids.nodes) * sizeof(Trie_Node));
    fprintf(stderr, "ids strings: #%zu, %zu\n", array_length(self->ids.strings),
        array_length(self->ids.strings) * sizeof(char));
    fprintf(stderr, "alleles: #%zu, %zu\n", array_length(self->alleles),
        array_length(self->alleles) * sizeof(Allele));
    fprintf(stderr, "join: #%zu, %zu\n", array_length(self->join), array_length(self->join) * sizeof(Join));

} // gva_index_stats
