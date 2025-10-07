#define _POSIX_C_SOURCE 200809L
#include <errno.h>      // errno
#include <limits.h>     // CHAR_BIT
#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // FILE, stderr, fclose, fopen, fprintf
#include <stdlib.h>     // EXIT_*

#include <string.h>     // strerror, strlen
#include <time.h>       // CLOCKS_PER_SEC, clock_t, clock

#include "../include/compare.h"     // bitset_fill
#include "../include/edit.h"        // gva_edit_distance
#include "../include/extractor.h"   // gva_canonical
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, GVA_Variant, gva_lcs_graph_*, gva_edges
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/string.h"      // GVA_String, gva_string_destroy
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "../include/utils.h"       // gva_fasta_sequence
#include "../include/variant.h"     // GVA_VARIANT_*, GVA_Variant, gva_parse_spdi
#include "array.h"          // ARRAY_DESTROY, array_length
#include "bitset.h"         // bitset_*
#include "common.h"         // MAX, MIN
#include "hash_table.h"     // HASH_TABLE_KEY, hash_table_*
#include "interval_tree.h"  // Interval_Tree, interval_tree_*
#include "trie.h"           // Trie, trie_*


#define TIC(x) clock_t const x = clock()
#define TOC(x) fprintf(stderr, "Elapsed time: %s (seconds): %.2f\n", (#x), (double)(clock() - (x)) / CLOCKS_PER_SEC)

#define LINE_SIZE 8194

#define REFERENCE_ID "NC_000001.11"


// FIXME: see `match_number` in src/variant.c
static inline size_t
parse_number(char const buffer[static 1], size_t idx[static 1])
{
    size_t number = 0;
    while (buffer[*idx] >= '0' && buffer[*idx] <= '9')
    {
        number = number * 10 + buffer[*idx] - '0';
        *idx += 1;
    } // while
    return number;
} // parse_number


GVA_Relation
compare_from_index(GVA_String const reference,
    GVA_Variant const lhs, gva_uint const lhs_dist,
    GVA_Variant const rhs, gva_uint const rhs_dist)
{
    // TODO: possibly reuse lhs_graph again

    GVA_Relation relation;

    if (gva_variant_eq(lhs, rhs))
    {
        return GVA_EQUIVALENT;
    } // if

    size_t const start = MIN(lhs.start, rhs.start);
    size_t const end = MAX(lhs.end, rhs.end);

    size_t const lhs_len = (lhs.start - start) + lhs.sequence.len + (end - lhs.end);
    size_t const rhs_len = (rhs.start - start) + rhs.sequence.len + (end - rhs.end);

    gva_uint distance = 0;
    if (lhs_len == 0)
    {
        distance = rhs_len;
    } // if
    else if (rhs_len == 0)
    {
        distance = lhs_len;
    } // if
    else
    {
        char* lhs_obs = gva_std_allocator.allocate(gva_std_allocator.context, NULL, 0, lhs_len);
        if (lhs_obs == NULL)
        {
            return -1;
        } // if
        memcpy(lhs_obs, reference.str + start, lhs.start - start);
        memcpy(lhs_obs + lhs.start - start, lhs.sequence.str, lhs.sequence.len);
        memcpy(lhs_obs + lhs.start - start + lhs.sequence.len, reference.str + lhs.end, end - lhs.end);

        char* rhs_obs = gva_std_allocator.allocate(gva_std_allocator.context, NULL, 0, rhs_len);
        if (rhs_obs == NULL)
        {
            lhs_obs = gva_std_allocator.allocate(gva_std_allocator.context, lhs_obs, lhs_len, 0);
            return -1;
        } // if
        memcpy(rhs_obs, reference.str + start, rhs.start - start);
        memcpy(rhs_obs + rhs.start - start, rhs.sequence.str, rhs.sequence.len);
        memcpy(rhs_obs + rhs.start - start + rhs.sequence.len, reference.str + rhs.end, end - rhs.end);

        distance = gva_edit_distance(gva_std_allocator, lhs_len, lhs_obs, rhs_len, rhs_obs);
        rhs_obs = gva_std_allocator.allocate(gva_std_allocator.context, rhs_obs, rhs_len, 0);
        lhs_obs = gva_std_allocator.allocate(gva_std_allocator.context, lhs_obs, lhs_len, 0);
    } // else

    if (lhs_dist + rhs_dist == distance)
    {
        return GVA_DISJOINT;
    } // if

    if (lhs_dist - rhs_dist == distance)
    {
        return GVA_CONTAINS;
    } // if

    if (rhs_dist - lhs_dist == distance)
    {
        return GVA_IS_CONTAINED;
    } // if

    size_t const len = end - start + 1;
    size_t const start_intersection = MAX(lhs.start, rhs.start);
    size_t const end_intersection = MIN(lhs.end, rhs.end);

    GVA_LCS_Graph lhs_graph = gva_lcs_graph_init(gva_std_allocator, lhs.end - lhs.start, reference.str + lhs.start, lhs.sequence.len, lhs.sequence.str, lhs.start);
    size_t* lhs_dels = bitset_init(gva_std_allocator, len);  // can be one shorter
    size_t* lhs_as = bitset_init(gva_std_allocator, len);
    size_t* lhs_cs = bitset_init(gva_std_allocator, len);
    size_t* lhs_gs = bitset_init(gva_std_allocator, len);
    size_t* lhs_ts = bitset_init(gva_std_allocator, len);
    bitset_fill(lhs_graph, start, start_intersection, end_intersection, lhs_dels, lhs_as, lhs_cs, lhs_gs, lhs_ts);

    GVA_LCS_Graph rhs_graph = gva_lcs_graph_init(gva_std_allocator, rhs.end - rhs.start, reference.str + rhs.start, rhs.sequence.len, rhs.sequence.str, rhs.start);
    size_t* rhs_dels = bitset_init(gva_std_allocator, len);  // can be one shorter
    size_t* rhs_as = bitset_init(gva_std_allocator, len);
    size_t* rhs_cs = bitset_init(gva_std_allocator, len);
    size_t* rhs_gs = bitset_init(gva_std_allocator, len);
    size_t* rhs_ts = bitset_init(gva_std_allocator, len);
    bitset_fill(rhs_graph, start, start_intersection, end_intersection, rhs_dels, rhs_as, rhs_cs, rhs_gs, rhs_ts);

    if (bitset_intersection_cnt(lhs_dels, rhs_dels) > 0 ||
        bitset_intersection_cnt(lhs_as, rhs_as) > 0 ||
        bitset_intersection_cnt(lhs_cs, rhs_cs) > 0 ||
        bitset_intersection_cnt(lhs_gs, rhs_gs) > 0 ||
        bitset_intersection_cnt(lhs_ts, rhs_ts) > 0)
    {
        relation = GVA_OVERLAP;
    } // if
    else
    {
        relation = GVA_DISJOINT;
    } // else

    lhs_ts = bitset_destroy(gva_std_allocator, lhs_ts);
    lhs_gs = bitset_destroy(gva_std_allocator, lhs_gs);
    lhs_cs = bitset_destroy(gva_std_allocator, lhs_cs);
    lhs_as = bitset_destroy(gva_std_allocator, lhs_as);
    lhs_dels = bitset_destroy(gva_std_allocator, lhs_dels);
    gva_lcs_graph_destroy(gva_std_allocator, lhs_graph);

    rhs_ts = bitset_destroy(gva_std_allocator, rhs_ts);
    rhs_gs = bitset_destroy(gva_std_allocator, rhs_gs);
    rhs_cs = bitset_destroy(gva_std_allocator, rhs_cs);
    rhs_as = bitset_destroy(gva_std_allocator, rhs_as);
    rhs_dels = bitset_destroy(gva_std_allocator, rhs_dels);
    gva_lcs_graph_destroy(gva_std_allocator, rhs_graph);

    return relation;
} // compare_from_index


static size_t
variants_distance(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref],
    GVA_Variant const lhs, GVA_Variant const rhs)
{
    size_t const start = MIN(lhs.start, rhs.start);
    size_t const end = MAX(lhs.end, rhs.end);

    size_t const lhs_len = (lhs.start - start) + lhs.sequence.len + (end - lhs.end);
    size_t const rhs_len = (rhs.start - start) + rhs.sequence.len + (end - rhs.end);

    if (lhs_len == 0)
    {
        return rhs_len;
    } // if
    if (rhs_len == 0)
    {
        return lhs_len;
    } // if

    char* lhs_obs = allocator.allocate(allocator.context, NULL, 0, lhs_len);
    char* rhs_obs = allocator.allocate(allocator.context, NULL, 0, rhs_len);
    if (lhs_obs == NULL || rhs_obs == NULL)
    {
        lhs_obs = allocator.allocate(allocator.context, lhs_obs, lhs_len, 0);
        rhs_obs = allocator.allocate(allocator.context, rhs_obs, rhs_len, 0);
        return -1;
    } // if

    memcpy(lhs_obs, reference + start, lhs.start - start);
    memcpy(lhs_obs + lhs.start - start, lhs.sequence.str, lhs.sequence.len);
    memcpy(lhs_obs + lhs.start - start + lhs.sequence.len, reference + lhs.end, end - lhs.end);

    memcpy(rhs_obs, reference + start, rhs.start - start);
    memcpy(rhs_obs + rhs.start - start, rhs.sequence.str, rhs.sequence.len);
    memcpy(rhs_obs + rhs.start - start + rhs.sequence.len, reference + rhs.end, end - rhs.end);

    size_t const distance = gva_edit_distance(allocator, lhs_len, lhs_obs, rhs_len, rhs_obs);

    rhs_obs = allocator.allocate(allocator.context, rhs_obs, rhs_len, 0);
    lhs_obs = allocator.allocate(allocator.context, lhs_obs, lhs_len, 0);

    return distance;
} // variants_distance


static GVA_Variant
construct_variant(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref],
    Interval_Tree const tree, Trie const trie, gva_uint* const nodes)
{
    size_t const n = array_length(nodes);
    GVA_Variant* variants = allocator.allocate(allocator.context, NULL, 0, n * sizeof(*variants));
    for (size_t i = 0; i < n; ++i)
    {
        variants[i] = (GVA_Variant) {
                tree.nodes[nodes[i]].start,
                tree.nodes[nodes[i]].end,
                trie_string(trie, tree.nodes[nodes[i]].inserted)
        };
    } // for

    GVA_Variant variant = {
            tree.nodes[nodes[0]].start,
            tree.nodes[nodes[n - 1]].end,
            {0, NULL}
    };
    variant.sequence = gva_string_concat(allocator, variant.sequence, variants[0].sequence);
    for (size_t i = 1; i < n; ++i)
    {
        variant.sequence = gva_string_concat(allocator, variant.sequence, (GVA_String) {variants[i].start - variants[i - 1].end, reference + variants[i - 1].end});
        variant.sequence = gva_string_concat(allocator, variant.sequence, variants[i].sequence);
    } // for
    // fprintf(stderr, GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(variant));

    // deallocate
    variants = allocator.allocate(allocator.context, variants, n * sizeof(*variants), 0);

    return variant;
} // construct_variant


static size_t
multiple_is_contained_distance(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref],
    GVA_LCS_Graph const graph,
    size_t const part_idx, Interval_Tree const tree, Trie const trie, gva_uint* const nodes)
{
    size_t const n = array_length(nodes);
    size_t lhs_distance = 0;
    for (size_t i = 0; i < n; ++i)
    {
        lhs_distance += tree.nodes[nodes[i]].distance;
    } // for

    size_t const rhs_distance = graph.local_supremal[part_idx + 1].distance;
    if (lhs_distance >= rhs_distance)
    {
        return 1;  // overlap
    } // if

    GVA_Variant const lhs = construct_variant(allocator, len_ref, reference, tree, trie, nodes);

    GVA_Variant rhs;
    gva_edges(graph.observed.str,
              graph.local_supremal[part_idx], graph.local_supremal[part_idx + 1],
              part_idx == 0, part_idx == array_length(graph.local_supremal) - 2,
              &rhs);

    size_t const distance = variants_distance(allocator, len_ref, reference, rhs, lhs);
    gva_string_destroy(allocator, lhs.sequence);

    if (rhs_distance - distance != lhs_distance)
    {
        return 1;  // overlap
    } // if
    return lhs_distance;
} // multiple_is_contained_distance


static inline size_t
prefix_length(size_t const len_lhs, char const lhs[static restrict len_lhs],
    size_t const len_rhs, char const rhs[static restrict len_rhs])
{
    size_t idx = 0;
    while (idx < len_lhs && idx < len_rhs && lhs[idx] == rhs[idx])
    {
        idx += 1;
    } // while
    return idx;
} // prefix_length


static inline size_t
suffix_length(size_t const len_lhs, char const lhs[static restrict len_lhs],
    size_t const len_rhs, char const rhs[static restrict len_rhs])
{
    size_t idx_lhs = len_lhs;
    size_t idx_rhs = len_rhs;
    while (idx_lhs > 0 && idx_rhs > 0 && lhs[idx_lhs - 1] == rhs[idx_rhs - 1])
    {
        idx_lhs -= 1;
        idx_rhs -= 1;
    } // while
    return len_lhs - idx_lhs;
} // suffix_length


static inline GVA_Variant
prefix_trimmed(size_t const len_ref, char const reference[static len_ref],
    GVA_Variant const variant)
{
    size_t const len = prefix_length((variant.end - variant.start), reference + variant.start, variant.sequence.len, variant.sequence.str);
    return (GVA_Variant) {variant.start + len, variant.end, {variant.sequence.len - len, variant.sequence.str + len}};
} // prefix_trimmed


static inline GVA_Variant
suffix_trimmed(size_t const len_ref, char const reference[static len_ref],
    GVA_Variant const variant)
{
    size_t const len = suffix_length((variant.end - variant.start), reference + variant.start, variant.sequence.len, variant.sequence.str);
    return (GVA_Variant) {variant.start, variant.end - len, {variant.sequence.len - len, variant.sequence.str}};
} // suffix_trimmed


int
vcf_main(int argc, char* argv[static argc + 1])
{
    errno = 0;
    FILE* stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    GVA_String reference = {0, NULL};
    reference = gva_fasta_sequence(gva_std_allocator, stream);
    fclose(stream);

    GVA_Variant variants[2];
    size_t count = 0;
    static char line[LINE_SIZE] = {0};
    while (fgets(line, sizeof(line), stdin) != NULL)
    {
        int const len = (char*) memchr(line, '\n', LINE_SIZE) - line;
        if (gva_parse_spdi(len, line, &variants[count > 0]) == 0)
        {
            fprintf(stderr, "error: SPDI parsing failed at line %zu: %s", count + 1, line);
            continue;
        } // if

        variants[count > 0] = prefix_trimmed(reference.len, reference.str, variants[count > 0]);

        for (size_t i = 0; i <= (count > 0); ++i)
        {
            fprintf(stderr, GVA_VARIANT_FMT " ", GVA_VARIANT_PRINT(variants[i]));
        } // for
        fprintf(stderr, "\n");

        GVA_LCS_Graph graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, 1 + (count > 0), variants);
        GVA_Variant local;
        for (size_t i = 0; i < array_length(graph.local_supremal) - 1; ++i)
        {
            gva_edges(graph.observed.str,
                graph.local_supremal[i], graph.local_supremal[i + 1],
                i == 0, i == array_length(graph.local_supremal) - 2,
                &local);
            if (i < array_length(graph.local_supremal) - 2)
            {
                fprintf(stderr, "OUT: " GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(local));
            } // if
        } // for
        gva_lcs_graph_destroy(gva_std_allocator, graph);

        if (count > 0)
        {
            gva_string_destroy(gva_std_allocator, variants[0].sequence);
        } // if
        variants[0] = (GVA_Variant) {local.start, local.end, gva_string_dup(gva_std_allocator, local.sequence)};
        variants[0] = suffix_trimmed(reference.len, reference.str, variants[0]);

        count += 1;
    } // while

    gva_string_destroy(gva_std_allocator, reference);
    return EXIT_SUCCESS;

    struct IN_EX
    {
        HASH_TABLE_KEY;
        gva_uint included;
        gva_uint excluded;
    }* table = hash_table_init(gva_std_allocator, 1024, sizeof(*table));

    HASH_TABLE_SET(gva_std_allocator, table, 42, ((struct IN_EX) {42, 1, 2}));

    if (table[HASH_TABLE_INDEX(table, 42)].gva_key == 42)
    {
        fprintf(stderr, "FOUND at %zu\n", HASH_TABLE_INDEX(table, 42));
    } // if

    table = HASH_TABLE_DESTROY(gva_std_allocator, table);

    return EXIT_SUCCESS;
}


int
main(int argc, char* argv[static argc + 1])
{
    errno = 0;
    FILE* stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    GVA_String reference = {0, NULL};
    reference = gva_fasta_sequence(gva_std_allocator, stream);
    fclose(stream);

    fprintf(stderr, "reference length: %zu\n", reference.len);

    Trie trie = trie_init();
    Interval_Tree tree = interval_tree_init();
    struct Node_Allele
    {
        gva_uint node;
        gva_uint allele;
        gva_uint next;
    }* node_allele_join = NULL;
    struct Allele
    {
        size_t   data;
        char*    spdi;  // FIXME
        gva_uint join_start;  // offset into node_allele_join
        gva_uint distance;
    }* db_alleles = NULL;

    size_t line_count = 0;
    static char line[LINE_SIZE] = {0};
    while (fgets(line, sizeof(line), stdin) != NULL)
    {
        size_t idx = 0;
        size_t const rsid = parse_number(line, &idx);
        idx += 1;  // skip space or tab
        int const len = (char*) memchr(line + idx, '\n', LINE_SIZE - idx) - (line + idx);
        GVA_Variant variant;
        if (gva_parse_spdi(len, line + idx, &variant) == 0)
        {
            fprintf(stderr, "error: SPDI parsing failed at line %zu: %s", line_count + 1, line);
            continue;
        } // if

        GVA_LCS_Graph const graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, 1, &variant);

        // FIXME
        size_t const spdi_len = snprintf(NULL, 0, GVA_VARIANT_FMT_SPDI, GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, graph.supremal)) + 1;
        char* const spdi = malloc(spdi_len);
        if (spdi == NULL)
        {
            fprintf(stderr, "spdi malloc() failed\n");
            return EXIT_FAILURE;
        }
        snprintf(spdi, spdi_len, GVA_VARIANT_FMT_SPDI, GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, graph.supremal));
        // end FIXME

        gva_uint const allele_idx = ARRAY_APPEND(
            gva_std_allocator, db_alleles,
            ((struct Allele) {rsid, spdi, array_length(node_allele_join), graph.distance})
        ) - 1;
        // fprintf(stderr, "input allele %zu: " GVA_VARIANT_FMT " (dist: %d) stored at: %d\n", rsid, GVA_VARIANT_PRINT(variant), graph.distance, allele_idx);

        for (size_t i = 0; i < array_length(graph.local_supremal) - 1; ++i)
        {
            GVA_Variant part;
            gva_edges(graph.observed.str,
                graph.local_supremal[i], graph.local_supremal[i + 1],
                i == 0, i == array_length(graph.local_supremal) - 2,
                &part);

            gva_uint const inserted_idx = trie_insert(gva_std_allocator, &trie, part.sequence.len, part.sequence.str);
            gva_uint const tmp_idx = ARRAY_APPEND(gva_std_allocator, tree.nodes, ((Interval_Tree_Node) {{GVA_NULL, GVA_NULL}, part.start, part.end, part.end, 0, inserted_idx, GVA_NULL, graph.local_supremal[i + 1].distance})) - 1;
            gva_uint const node_idx = interval_tree_insert(&tree, tmp_idx);
            if (node_idx != tmp_idx)
            {
                // fprintf(stderr, "    part %zu already in the tree\n", i);
                array_header(tree.nodes)->length -= 1;  // reverse; node already in the tree
            } // if
            tree.nodes[node_idx].alleles = ARRAY_APPEND(gva_std_allocator, node_allele_join, ((struct Node_Allele) {node_idx, allele_idx, tree.nodes[node_idx].alleles})) - 1;
            // fprintf(stderr, "    input part %zu (allele %d): " GVA_VARIANT_FMT " stored at: %d\n", i, allele_idx, GVA_VARIANT_PRINT(part), node_idx);
        } // for

        gva_string_destroy(gva_std_allocator, graph.observed);
        gva_lcs_graph_destroy(gva_std_allocator, graph);

        line_count += 1;
    } // while
    fprintf(stderr, "line count: %zu\n", line_count);

    fprintf(stderr, "tree nodes: %zu\n", array_length(tree.nodes));
    fprintf(stderr, "trie nodes: %zu (%zu)\n", array_length(trie.nodes), trie.strings.len);

    fprintf(stderr, "#db_alleles: %zu\n", array_length(db_alleles));
    fprintf(stderr, "#join:    %zu\n", array_length(node_allele_join));

    if (db_alleles == NULL || node_allele_join == NULL)
    {
        return -1;
    } // if

    errno = 0;
    stream = fopen(argv[2], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    // for every query
    line_count = 0;
    while (fgets(line, sizeof(line), stream) != NULL)
    {
        size_t idx = 0;
        size_t const query_id = parse_number(line, &idx);
        idx += 1;  // skip space or tab
        int const len = (char*) memchr(line + idx, '\n', LINE_SIZE - idx) - (line + idx);
        GVA_Variant query_var;
        if (gva_parse_spdi(len, line + idx, &query_var) == 0)
        {
            fprintf(stderr, "error: SPDI parsing failed at line %zu: %s", line_count + 1, line);
            continue;
        } // if

        // Join nodes in the index to parts in the query
        struct NODE_PARTS
        {
            HASH_TABLE_KEY;
            GVA_Relation relation;
            gva_uint start;
            gva_uint end;
            gva_uint included;
        }* node_parts_table = hash_table_init(gva_std_allocator, 1024, sizeof(*node_parts_table));

        // query graph
        GVA_LCS_Graph const query_graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, 1, &query_var);
        // fprintf(stderr, "query allele %zu: " GVA_VARIANT_FMT " (dist: %u)\n", query_id, GVA_VARIANT_PRINT(query_var), query_graph.distance);

        for (size_t part_idx = 0; part_idx < array_length(query_graph.local_supremal) - 1; ++part_idx)
        {
            // query local supremal part
            GVA_Variant query_part;
            gva_edges(query_graph.observed.str,
                      query_graph.local_supremal[part_idx], query_graph.local_supremal[part_idx + 1],
                      part_idx == 0, part_idx == array_length(query_graph.local_supremal) - 2,
                      &query_part);

            gva_uint const query_dist = query_graph.local_supremal[part_idx + 1].distance;
            // fprintf(stderr, "    query part %zu (dist: %d): " GVA_VARIANT_FMT "\n", part_idx, query_dist, GVA_VARIANT_PRINT(query_part));

            // find candidates for this local supremal part
            gva_uint* candidates = interval_tree_intersection(gva_std_allocator, tree, query_part.start, query_part.end);
            for (size_t can_idx = 0; can_idx < array_length(candidates); ++can_idx)
            {
                gva_uint const node_idx = candidates[can_idx];
                GVA_Variant const db_var = {tree.nodes[node_idx].start, tree.nodes[node_idx].end,
                                            trie_string(trie, tree.nodes[node_idx].inserted)};
                GVA_Relation const relation = compare_from_index(reference, db_var, tree.nodes[node_idx].distance, query_part, query_dist);
                // fprintf(stderr, "        db candidate node_idx: %d (dist: %d): " GVA_VARIANT_FMT " %s\n", node_idx, tree.nodes[node_idx].distance, GVA_VARIANT_PRINT(db_var), GVA_RELATION_LABELS[relation]);

                if (relation == GVA_DISJOINT)
                {
                    continue;
                }

                // link nodes to query parts
                size_t hash_idx = HASH_TABLE_INDEX(node_parts_table, node_idx);
                if (node_parts_table[hash_idx].gva_key != node_idx)
                {
                    HASH_TABLE_SET(gva_std_allocator, node_parts_table, node_idx,
                                   ((struct NODE_PARTS) {node_idx, relation, part_idx, part_idx + 1, query_dist}));
                    hash_idx = HASH_TABLE_INDEX(node_parts_table, node_idx);
                } // if

                if (relation == GVA_EQUIVALENT || relation == GVA_IS_CONTAINED)
                {
                    node_parts_table[hash_idx].included = tree.nodes[node_idx].distance;
                } // if
                else if (relation == GVA_CONTAINS)
                {
                    node_parts_table[hash_idx].end = part_idx + 1;
                } // if
                else // GVA_OVERLAP
                {
                    node_parts_table[hash_idx].included = 1;
                } // else
            } // for every candidate
            candidates = ARRAY_DESTROY(gva_std_allocator, candidates);
        } // for query allele parts

        // We have now compared all query parts to the index

        // fix containment for multiple parts in single node for every query
        for (size_t npt_index = 0; npt_index < array_header(node_parts_table)->capacity; ++npt_index)
        {
            size_t const node_idx = node_parts_table[npt_index].gva_key;
            if (node_idx == (uint32_t) -1)
            {
                continue;
            } // if

            if (node_parts_table[npt_index].relation == GVA_CONTAINS &&
                node_parts_table[npt_index].end - node_parts_table[npt_index].start > 1)
            {
                // fprintf(stderr, "    Perform extra comparison\n");
                size_t slice_dist = 0;
                for (size_t i = node_parts_table[npt_index].start; i < node_parts_table[npt_index].end; ++i)
                {
                    sum_distance += query_graph.local_supremal[i + 1].distance;
                } // for
                // fprintf(stderr, "    slice dist: %zu node dist: %u\n", slice_dist, tree.nodes[node_idx].distance);
                if (slice_dist >= tree.nodes[node_idx].distance)
                {
                    // fprintf(stderr, "here A\n");
                    node_parts_table[npt_index].relation = GVA_OVERLAP;
                    node_parts_table[npt_index].included = 1;
                    continue;
                } // if

                GVA_Variant const lhs = {tree.nodes[node_idx].start,
                                         tree.nodes[node_idx].end,
                                         trie_string(trie, tree.nodes[node_idx].inserted)};

                GVA_Variant rhs;
                gva_edges(query_graph.observed.str,
                          query_graph.local_supremal[node_parts_table[npt_index].start],
                          query_graph.local_supremal[node_parts_table[npt_index].end],
                          node_parts_table[npt_index].start == 0,
                          node_parts_table[npt_index].end == array_length(query_graph.local_supremal) - 1,
                          &rhs);
                // fprintf(stderr, GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(rhs));

                size_t const op_distance = variants_distance(gva_std_allocator, reference.len, reference.str, lhs, rhs);
                // fprintf(stderr, "op_distance: %d\n", op_distance);

                node_parts_table[npt_index].included = slice_dist;

                if (tree.nodes[node_idx].distance - op_distance != slice_dist)
                {
                    // fprintf(stderr, "here B\n");
                    node_parts_table[npt_index].relation = GVA_OVERLAP;
                    node_parts_table[npt_index].included = 1;
                    // fprintf(stderr, "OVERLAP!!\n");
                } // if
            } // if
        } // for node_parts_table

        // fprintf(stderr, "db_allele table:\n");
        // for (size_t i = 0; i < array_length(db_alleles); ++i)
        // {
        //     fprintf(stderr, "entry: %zu %s %zu %zu\n",
        //             db_alleles[i].data,
        //             db_alleles[i].spdi,
        //             db_alleles[i].join_start,
        //             db_alleles[i].distance
        //     );
        // }
        // fprintf(stderr, "\n");
        //
        // fprintf(stderr, "node_allele_join table:\n");
        // for (size_t i = 0; i < array_length(node_allele_join); ++i)
        // {
        //     fprintf(stderr, "entry: %zu %zu %d\n",
        //             node_allele_join[i].node,
        //             node_allele_join[i].allele,
        //             node_allele_join[i].next
        //     );
        // }
        // fprintf(stderr, "\n");

        struct RESULT_ALLELES
        {
            HASH_TABLE_KEY;
        }* results_table = hash_table_init(gva_std_allocator, 1024, sizeof(*results_table));

        // fprintf(stderr, "Loop over all node parts for every query:\n");
        for (size_t npt_index = 0; npt_index < array_header(node_parts_table)->capacity; ++npt_index)
        {
            size_t const node_idx = node_parts_table[npt_index].gva_key;
            if (node_idx == (uint32_t) - 1)
            {
                continue;
            } // if

            for (size_t naj_table_idx = tree.nodes[node_idx].alleles;
                 naj_table_idx != GVA_NULL;
                 naj_table_idx = node_allele_join[naj_table_idx].next)
            {
                size_t const allele_idx = node_allele_join[naj_table_idx].allele;
                // TODO: direct SET
                size_t hash_idx = HASH_TABLE_INDEX(results_table, allele_idx);
                if (results_table[hash_idx].gva_key != allele_idx)
                {
                    HASH_TABLE_SET(gva_std_allocator, results_table, allele_idx, ((struct RESULT_ALLELES) {allele_idx}));
                } // if
                // TODO: construct allele distance here?!
            } // for alleles
        } // for node_parts_table

        // build result vector for every query
        for (size_t results_index = 0; results_index < array_header(results_table)->capacity; ++results_index)
        {
            size_t allele_idx = results_table[results_index].gva_key;
            if (allele_idx == (uint32_t) - 1)
            {
                continue;
            } // if
            // fprintf(stderr, "allele_idx: %zu\n", allele_idx);

            GVA_Relation relation = GVA_DISJOINT;
            gva_uint included = 0;

            gva_uint* is_contained_nodes = NULL;
            gva_uint is_contained_part_idx = -1;

            // loop over all nodes for this allele
            // fprintf(stderr, "Loop over nodes for this allele:\n");
            for (size_t join_idx = db_alleles[allele_idx].join_start;
                join_idx < array_length(node_allele_join) && node_allele_join[join_idx].allele == allele_idx;
                ++join_idx)
            {
                // fprintf(stderr, "join_idx: %zu\n", join_idx);
                size_t const node_idx = node_allele_join[join_idx].node;
                size_t const hash_idx = HASH_TABLE_INDEX(node_parts_table, node_idx);
                if (node_idx != node_parts_table[hash_idx].gva_key)
                {
                    // fprintf(stderr, "when?!?\n");
                    continue;
                } // if
                // fprintf(stderr, "node_idx: %zu\n", node_idx);

                if (node_parts_table[hash_idx].relation == GVA_EQUIVALENT)
                {
                    included += node_parts_table[hash_idx].included;

                    if (relation == GVA_EQUIVALENT || relation == GVA_DISJOINT)
                    {
                        relation = GVA_EQUIVALENT;
                    }
                } // if
                else if (node_parts_table[hash_idx].relation == GVA_CONTAINS)
                {
                    // fprintf(stderr, "here D\n");
                    included += node_parts_table[hash_idx].included;

                    if (relation == GVA_IS_CONTAINED)
                    {
                        // included = 1; can it be an overestimate ?!?
                        relation = GVA_OVERLAP;
                        break;
                    }
                    relation = GVA_CONTAINS;
                }
                else if (node_parts_table[hash_idx].relation == GVA_IS_CONTAINED)
                {
                    size_t const part_idx = node_parts_table[hash_idx].start;
                    // fprintf(stderr, "part_idx: %zu\n", part_idx);
                    // fprintf(stderr, "if containment current: %s\n", GVA_RELATION_LABELS[relation]);
                    if (relation == GVA_CONTAINS)
                    {
                        included += node_parts_table[hash_idx].included;  // TODO: included = 1; ???
                        relation = GVA_OVERLAP;
                        break;
                    } // if

                    if (part_idx != is_contained_part_idx)
                    {
                        // fprintf(stderr, "close window upper %s\n", GVA_RELATION_LABELS[relation]);
                        // fprintf(stderr, "#is_contained_nodes: %zu\n", array_length(is_contained_nodes));
                        // close old window
                        if (array_length(is_contained_nodes) > 1)
                        {
                            // fprintf(stderr, "Found multiple is_containeds for the same part\n");
                            size_t const slice_dist = multiple_is_contained_distance(gva_std_allocator, reference.len, reference.str,
                                                                                     query_graph, part_idx, tree, trie, is_contained_nodes);
                            if (slice_dist == 1)
                            {
                                included = 1;
                                relation = GVA_OVERLAP;
                                break;
                            } // if
                            included += slice_dist;
                        } // if
                        else if (array_length(is_contained_nodes) == 1)
                        {
                            // fprintf(stderr, "Now should add previous node\n");
                            included += tree.nodes[is_contained_nodes[0]].distance;
                        }

                        // open new window
                        // fprintf(stderr, "open new window\n");
                        is_contained_part_idx = part_idx;
                        is_contained_nodes = ARRAY_DESTROY(gva_std_allocator, is_contained_nodes);
                    } // if part_idx != is_contained_part_idx
                    ARRAY_APPEND(gva_std_allocator, is_contained_nodes, node_idx);
                    relation = GVA_IS_CONTAINED;
                } // if relation == GVA_IS_CONTAINED
                else if (node_parts_table[hash_idx].relation == GVA_OVERLAP)
                {
                    // fprintf(stderr, "this overlap is reached\n");
                    included = 1;
                    relation = GVA_OVERLAP;
                    break;
                } // if
                // fprintf(stderr, "here F\n");

            } // for all nodes for this allele

            // fprintf(stderr, "here G\n");

            // fprintf(stderr, "included: %zu\n", included);

            // close old window
            size_t const n = array_length(is_contained_nodes);
            if (relation == GVA_IS_CONTAINED && n == 1)
            {
                // fprintf(stderr, "found single is_contained\n");
                included += tree.nodes[is_contained_nodes[0]].distance;
                // fprintf(stderr, "increase included with %u to %zu\n", tree.nodes[is_contained_nodes[0]].distance, included);
            }
            else if (relation == GVA_IS_CONTAINED && n > 1)
            {
                // fprintf(stderr, "Found multiple is_containeds for the same part at the end\n");
                size_t const slice_dist = multiple_is_contained_distance(gva_std_allocator, reference.len, reference.str,
                                                                         query_graph, is_contained_part_idx, tree, trie, is_contained_nodes);
                if (slice_dist == 1)
                {
                    included = 1;
                    relation = GVA_OVERLAP;
                } // if
                else
                {
                    included += slice_dist;
                    // fprintf(stderr, "increase included with %zu to %zu\n", slice_dist, included);
                } // else

            } // if is_contained repair
            is_contained_nodes = ARRAY_DESTROY(gva_std_allocator, is_contained_nodes);

            if (included > 0)
            {
                // fprintf(stderr, "allele dist: %u query dist: %u\n", db_alleles[allele_idx].distance, query_graph.distance);
                // fprintf(stderr, "relation: %s\n", GVA_RELATION_LABELS[relation]);

                gva_uint const allele_excluded = db_alleles[allele_idx].distance - included;
                gva_uint const query_excluded = query_graph.distance - included;

                // fprintf(stderr, "included: %zu, allele_excluded: %d query_excluded: %zu\n",
                //         included, allele_excluded, query_excluded);
                if (allele_excluded > 0 && query_excluded > 0)
                {
                    relation = GVA_OVERLAP;
                }
                else if (allele_excluded > 0)
                {
                    relation = GVA_CONTAINS;
                }
                else if (query_excluded > 0)
                {
                    relation = GVA_IS_CONTAINED;
                }
                else // if (allele_excluded == 0 && query_excluded == 0)
                {
                    relation = GVA_EQUIVALENT;
                }
                // fprintf(stderr, "allele_idx: %zu relation: %s in: %zu\n", allele_idx, GVA_RELATION_LABELS[relation], included);

                // only for testing
                if (relation != GVA_EQUIVALENT || db_alleles[allele_idx].data < query_id)
                {
                    printf("%zu %s %zu " GVA_VARIANT_FMT_SPDI " %s\n",
                            db_alleles[allele_idx].data,
                            db_alleles[allele_idx].spdi,
                            query_id,
                            GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, query_graph.supremal),
                            GVA_RELATION_LABELS[relation]);
                }
            } // if
        } // for all alleles


        // fprintf(stderr, "\n");

        results_table = HASH_TABLE_DESTROY(gva_std_allocator, results_table);
        node_parts_table = HASH_TABLE_DESTROY(gva_std_allocator, node_parts_table);

        gva_string_destroy(gva_std_allocator, query_graph.observed);
        gva_lcs_graph_destroy(gva_std_allocator, query_graph);

        // fprintf(stderr, "\n");
        line_count += 1;
    } // while fgets
    fclose(stream);

    for (size_t allele_idx = 0; allele_idx < array_length(db_alleles); ++allele_idx)
    {
        free(db_alleles[allele_idx].spdi);
    }
    db_alleles = ARRAY_DESTROY(gva_std_allocator, db_alleles);
    node_allele_join = ARRAY_DESTROY(gva_std_allocator, node_allele_join);
    interval_tree_destroy(gva_std_allocator, &tree);
    trie_destroy(gva_std_allocator, &trie);

    gva_string_destroy(gva_std_allocator, reference);

    return EXIT_SUCCESS;
} // main
