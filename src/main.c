#define _POSIX_C_SOURCE 200809L
#include <errno.h>      // errno
#include <limits.h>     // CHAR_BIT
#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // FILE, stderr, stdout, fclose, fopen, fprintf
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
#include "hash_table.h"     // GVA_NOT_FOUND, HASH_TABLE_KEY, hash_table_*
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


GVA_String
gva_fasta_sequence_blob(GVA_Allocator const allocator, FILE* const stream)
{
    GVA_String reference = {0, NULL};

    static char line[42] = {0};
    if (fgets(line, sizeof(line), stream) == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return reference;
    } // if
    size_t idx = 0;
    size_t length = parse_number(line, &idx);

    reference.str = allocator.allocate(allocator.context, NULL, 0, length);
    if (reference.str == NULL)
    {
        fprintf(stderr, "allocate error: %s\n", strerror(errno));
        return reference;
    } // if

    if (fgets((char *)reference.str, length, stream) == NULL)
    {
        fprintf(stderr, "read error: %s\n", strerror(errno));
        gva_string_destroy(gva_std_allocator, reference);
        return reference;
    } // if
    reference.len = length;

    return reference;
} // gva_fasta_sequence_blob


int
fasta_blob_write(int argc, char* argv[static argc + 1])
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

    fprintf(stdout, "%zu\n", reference.len);
    size_t const ret = fwrite(reference.str, 1, reference.len, stdout);
    if (ret != reference.len)
    {
        return EXIT_FAILURE;
    } // if

    return EXIT_SUCCESS;
} // fasta_blob_write


int
fasta_blob_read(int argc, char* argv[static argc + 1])
{
    errno = 0;
    FILE* stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    GVA_String reference = gva_fasta_sequence_blob(gva_std_allocator, stream);
    fclose(stream);

    fprintf(stderr, "reference length: %zu\n", reference.len);
    gva_string_destroy(gva_std_allocator, reference);

    return EXIT_SUCCESS;
} // fasta_blob_read


GVA_Relation
compare_from_index(GVA_String const reference,
    GVA_Variant const lhs, gva_uint const lhs_distance,
    GVA_Variant const rhs, gva_uint const rhs_distance)
{
    // TODO: possibly reuse rhs_graph

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

    if (lhs_distance + rhs_distance == distance)
    {
        return GVA_DISJOINT;
    } // if

    if (lhs_distance - rhs_distance == distance)
    {
        return GVA_CONTAINS;
    } // if

    if (rhs_distance - lhs_distance == distance)
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

    GVA_Relation relation = GVA_DISJOINT;
    if (bitset_intersection_cnt(lhs_dels, rhs_dels) > 0 ||
        bitset_intersection_cnt(lhs_as, rhs_as) > 0 ||
        bitset_intersection_cnt(lhs_cs, rhs_cs) > 0 ||
        bitset_intersection_cnt(lhs_gs, rhs_gs) > 0 ||
        bitset_intersection_cnt(lhs_ts, rhs_ts) > 0)
    {
        relation = GVA_OVERLAP;
    } // if

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
    // FIXME: in 1 loop
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


int vcf_main(int argc, char* argv[static argc + 1])
{
    if (argc < 3)
    {
        fprintf(stderr, "usage: %s reference.blob gap\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    errno = 0;
    FILE *stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    GVA_String reference = gva_fasta_sequence_blob(gva_std_allocator, stream);
    fclose(stream);
    fprintf(stderr, "reference length: %zu\n", reference.len);

    size_t const gap = atoll(argv[2]);
    fprintf(stderr, "gap: %zu\n", gap);

    GVA_Variant* variants = NULL;
    GVA_Variant* lss = NULL;

    size_t last = 0;
    size_t line_count = 0;
    size_t dropped = 0;
    size_t distance = 0;
    static char line[LINE_SIZE] = {0};
    while (fgets(line, sizeof(line), stdin) != NULL)
    {
        //size_t idx = 0;
        //size_t const distance = parse_number(line, &idx);
        //idx += 1;  // skip space or tab
        //int const len = (char*) memchr(line + idx, '\n', LINE_SIZE - idx) - (line + idx);
        size_t const len = strlen(line) - 1;
        GVA_Variant variant;
        if (gva_parse_spdi(len, line, &variant) == 0)
        {
            fprintf(stderr, "error: SPDI parsing failed at line %zu: %s", line_count + 1, line);
            continue;
        } // if

        GVA_Variant trimmed = prefix_trimmed(reference.len, reference.str, variant);

        if (array_length(variants) > 0 && trimmed.start < variants[array_length(variants) - 1].end)
        {
            //fprintf(stderr, "dropped: at line %zu: %s", line_count + 1, line);
            dropped += 1;
            continue;
        } // if

        if (array_length(variants) > 0 && trimmed.start - variants[array_length(variants) - 1].end > gap)
        {
            //fprintf(stderr, "Maak graaf: %zu\n", last);
            GVA_LCS_Graph graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, array_length(variants) - last, variants + last);

            for (size_t i = 0; i < array_length(graph.local_supremal) - 1; ++i)
            {
                GVA_Variant variant;
                gva_edges(graph.observed.str,
                      graph.local_supremal[i], graph.local_supremal[i + 1],
                      i == 0, i == array_length(graph.local_supremal) - 2, &variant);
                fprintf(stdout, "%u " GVA_VARIANT_FMT_SPDI "\n", graph.local_supremal[i + 1].distance, GVA_VARIANT_PRINT_SPDI("NC_000001.11", variant));
                distance += graph.local_supremal[i + 1].distance;
                ARRAY_APPEND(gva_std_allocator, lss, gva_variant_dup(gva_std_allocator, variant));
            } // for

            last = array_length(variants);

            gva_string_destroy(gva_std_allocator, graph.observed);
            gva_lcs_graph_destroy(gva_std_allocator, graph);
        } // if

        //fprintf(stderr, "Add variant\n");

        ARRAY_APPEND(gva_std_allocator, variants, gva_variant_dup(gva_std_allocator, trimmed));
        line_count += 1;
    } // while
    if (last < array_length(variants))
    {
        //fprintf(stderr, "Maak graaf: %zu\n", last);
        GVA_LCS_Graph graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, array_length(variants) - last, variants + last);

        for (size_t i = 0; i < array_length(graph.local_supremal) - 1; ++i)
        {
            GVA_Variant variant;
            gva_edges(graph.observed.str,
                  graph.local_supremal[i], graph.local_supremal[i + 1],
                  i == 0, i == array_length(graph.local_supremal) - 2, &variant);
            fprintf(stdout, "%u " GVA_VARIANT_FMT_SPDI "\n", graph.local_supremal[i + 1].distance, GVA_VARIANT_PRINT_SPDI("NC_000001.11", variant));
            distance += graph.local_supremal[i + 1].distance;
            ARRAY_APPEND(gva_std_allocator, lss, gva_variant_dup(gva_std_allocator, variant));
        } // for

        last = array_length(variants);

        gva_string_destroy(gva_std_allocator, graph.observed);
        gva_lcs_graph_destroy(gva_std_allocator, graph);
    } // if

    GVA_String observed = gva_patch(gva_std_allocator, reference.len, reference.str, array_length(variants), variants);
    GVA_String lss_observed = gva_patch(gva_std_allocator, reference.len, reference.str, array_length(lss), lss);
    size_t const in_distance = gva_edit_distance(gva_std_allocator, reference.len, reference.str, observed.len, observed.str);

    fprintf(stderr, "===INPUT===\n");
    fprintf(stderr, "#variants: %zu\n", array_length(variants));
    fprintf(stderr, "#dropped:  %zu\n", dropped);
    fprintf(stderr, "distance:  %zu\n", in_distance);
    fprintf(stderr, "===OUTPUT===\n");
    fprintf(stderr, "#parts:    %zu\n", array_length(lss));
    fprintf(stderr, "distance:  %zu\n", distance);

    for (size_t i = 0; i < array_length(variants); ++i)
    {
        gva_string_destroy(gva_std_allocator, variants[i].sequence);
    } // for
    variants = ARRAY_DESTROY(gva_std_allocator, variants);

    for (size_t i = 0; i < array_length(lss); ++i)
    {
        gva_string_destroy(gva_std_allocator, lss[i].sequence);
    } // for
    lss = ARRAY_DESTROY(gva_std_allocator, lss);

    bool const eq = observed.len == lss_observed.len && strncmp(observed.str, lss_observed.str, observed.len) == 0;

    gva_string_destroy(gva_std_allocator, lss_observed);
    gva_string_destroy(gva_std_allocator, observed);
    gva_string_destroy(gva_std_allocator, reference);

    if (eq)
    {
        return EXIT_SUCCESS;
    } // if

    fprintf(stderr, "ERROR\n");

    return EXIT_FAILURE;
} // vcf_main


void
repair_is_contained(GVA_Allocator const allocator, GVA_String const reference, GVA_LCS_Graph const graph, Interval_Tree const tree, Trie const trie,
      gva_uint** nodes, gva_uint const part_idx, gva_uint* const included, GVA_Relation* const relation)
{
    size_t const nodes_len = array_length(*nodes);
    if (*relation == GVA_IS_CONTAINED && nodes_len == 1)
    {
        *included += tree.nodes[*nodes[0]].distance;
    } // if
    else if (*relation == GVA_IS_CONTAINED && nodes_len > 1)
    {
        size_t const slice_dist = multiple_is_contained_distance(allocator, reference.len, reference.str,
                                                                 graph, part_idx, tree, trie, *nodes);
        if (slice_dist == 1)
        {
            *included = 1;
            *relation = GVA_OVERLAP;
        } // if
        else
        {
            *included += slice_dist;
        } // else
    } // if
    *nodes = ARRAY_DESTROY(gva_std_allocator, *nodes);
} // repair_is_contained


int
dbsnp_main(int argc, char* argv[static argc + 1])
{
    errno = 0;
    FILE* stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    GVA_String reference = {0, NULL};
    reference = gva_fasta_sequence_blob(gva_std_allocator, stream);
    fclose(stream);

    fprintf(stderr, "reference length: %zu\n", reference.len);

    Trie trie = trie_init();
    Interval_Tree tree = interval_tree_init();
    struct Node_Allele
    {
        gva_uint link;
        gva_uint next;
    }* node_allele_join = NULL;
    struct Allele
    {
        gva_uint line;
        gva_uint join_start;  // offset into node_allele_join
        gva_uint join_end;  // end of offset into node_allele_join
        gva_uint distance;
    }* db_alleles = NULL;

    size_t line_count = 0;
    static char line[LINE_SIZE] = {0};
    while (fgets(line, sizeof(line), stdin) != NULL)
    {
        size_t idx = 0;
        parse_number(line, &idx);
        idx += 1;  // skip space or tab
        int const len = (char*) memchr(line + idx, '\n', LINE_SIZE - idx) - (line + idx);
        GVA_Variant variant;
        if (gva_parse_spdi(len, line + idx, &variant) == 0)
        {
            fprintf(stderr, "error: SPDI parsing failed at line %zu: %s", line_count + 1, line);
            continue;
        } // if

        GVA_LCS_Graph const graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, 1, &variant);

        gva_uint const allele_idx = ARRAY_APPEND(
            gva_std_allocator, db_alleles,
            ((struct Allele) {line_count, array_length(node_allele_join),
                              array_length(node_allele_join) + array_length(graph.local_supremal) - 1,
                    graph.distance})
        ) - 1;

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
                array_header(tree.nodes)->length -= 1;  // reverse; node already in the tree
            } // if
            tree.nodes[node_idx].alleles = ARRAY_APPEND(gva_std_allocator, node_allele_join, ((struct Node_Allele) {node_idx ^ allele_idx, tree.nodes[node_idx].alleles})) - 1;
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
        parse_number(line, &idx);
        idx += 1;  // skip space or tab
        int const len = (char*) memchr(line + idx, '\n', LINE_SIZE - idx) - (line + idx);
        GVA_Variant rhs_var;
        if (gva_parse_spdi(len, line + idx, &rhs_var) == 0)
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

        GVA_LCS_Graph const rhs_graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, 1, &rhs_var);

        for (size_t part_idx = 0; part_idx < array_length(rhs_graph.local_supremal) - 1; ++part_idx)
        {
            GVA_Variant rhs_part;
            gva_edges(rhs_graph.observed.str,
                      rhs_graph.local_supremal[part_idx], rhs_graph.local_supremal[part_idx + 1],
                      part_idx == 0, part_idx == array_length(rhs_graph.local_supremal) - 2,
                      &rhs_part);

            gva_uint const rhs_distance = rhs_graph.local_supremal[part_idx + 1].distance;

            gva_uint* candidates = interval_tree_intersection(gva_std_allocator, tree, rhs_part.start, rhs_part.end);
            for (size_t can_idx = 0; can_idx < array_length(candidates); ++can_idx)
            {
                gva_uint const node_idx = candidates[can_idx];
                GVA_Variant const db_var = {tree.nodes[node_idx].start, tree.nodes[node_idx].end,
                                            trie_string(trie, tree.nodes[node_idx].inserted)};
                GVA_Relation const relation = compare_from_index(reference, db_var, tree.nodes[node_idx].distance, rhs_part, rhs_distance);
                if (relation == GVA_DISJOINT)
                {
                    continue;
                } // if

                // link nodes to query parts
                size_t hash_idx = HASH_TABLE_INDEX(node_parts_table, node_idx);
                if (node_parts_table[hash_idx].gva_key != node_idx)
                {
                    HASH_TABLE_SET(gva_std_allocator, node_parts_table, node_idx,
                                   ((struct NODE_PARTS) {node_idx, relation, part_idx, part_idx + 1, rhs_distance}));
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
            if (node_idx == GVA_NOT_FOUND)
            {
                continue;
            } // if

            if (node_parts_table[npt_index].relation == GVA_CONTAINS &&
                node_parts_table[npt_index].end - node_parts_table[npt_index].start > 1)
            {
                size_t const lhs_distance = tree.nodes[node_idx].distance;
                size_t rhs_distance = 0;
                for (size_t i = node_parts_table[npt_index].start; i < node_parts_table[npt_index].end; ++i)
                {
                    rhs_distance += rhs_graph.local_supremal[i + 1].distance;
                } // for
                if (rhs_distance >= lhs_distance)
                {
                    node_parts_table[npt_index].included = 1;
                    node_parts_table[npt_index].relation = GVA_OVERLAP;
                    continue;
                } // if

                GVA_Variant const lhs = {tree.nodes[node_idx].start,
                                         tree.nodes[node_idx].end,
                                         trie_string(trie, tree.nodes[node_idx].inserted)};

                GVA_Variant rhs;
                gva_edges(rhs_graph.observed.str,
                          rhs_graph.local_supremal[node_parts_table[npt_index].start],
                          rhs_graph.local_supremal[node_parts_table[npt_index].end],
                          node_parts_table[npt_index].start == 0,
                          node_parts_table[npt_index].end == array_length(rhs_graph.local_supremal) - 1,
                          &rhs);

                size_t const distance = variants_distance(gva_std_allocator, reference.len, reference.str, lhs, rhs);
                if (lhs_distance - distance == rhs_distance)
                {
                    node_parts_table[npt_index].included = rhs_distance;
                } // if
                else
                {
                    node_parts_table[npt_index].included = 1;
                    node_parts_table[npt_index].relation = GVA_OVERLAP;
                } // if
            } // if
        } // for node_parts_table

        struct RESULT_ALLELES
        {
            HASH_TABLE_KEY;
        }* results_table = hash_table_init(gva_std_allocator, 1024, sizeof(*results_table));

        // find all alleles that were part of a non-disjoint relation
        for (size_t npt_idx = 0; npt_idx < array_header(node_parts_table)->capacity; ++npt_idx)
        {
            size_t const node_idx = node_parts_table[npt_idx].gva_key;
            if (node_idx == GVA_NOT_FOUND)
            {
                continue;
            } // if

            for (size_t naj_table_idx = tree.nodes[node_idx].alleles;
                 naj_table_idx != GVA_NULL;
                 naj_table_idx = node_allele_join[naj_table_idx].next)
            {
                size_t const allele_idx = node_allele_join[naj_table_idx].link ^ node_idx;
                HASH_TABLE_SET(gva_std_allocator, results_table, allele_idx, ((struct RESULT_ALLELES) {allele_idx}));
            } // for alleles
        } // for node_parts_table

        // build result vector for every query
        for (size_t results_idx = 0; results_idx < array_header(results_table)->capacity; ++results_idx)
        {
            size_t allele_idx = results_table[results_idx].gva_key;
            if (allele_idx == GVA_NOT_FOUND)
            {
                continue;
            } // if

            gva_uint included = 0;
            GVA_Relation relation = GVA_DISJOINT;

            gva_uint* is_contained_nodes = NULL;
            gva_uint is_contained_part_idx = -1;

            // loop over all nodes for this allele
            for (size_t join_idx = db_alleles[allele_idx].join_start; join_idx < db_alleles[allele_idx].join_end; ++join_idx)
            {
                size_t node_idx = node_allele_join[join_idx].link ^ allele_idx;
                size_t const hash_idx = HASH_TABLE_INDEX(node_parts_table, node_idx);
                if (node_idx != node_parts_table[hash_idx].gva_key)
                {
                    continue;
                } // if

                if (node_parts_table[hash_idx].relation == GVA_EQUIVALENT)
                {
                    included += node_parts_table[hash_idx].included;
                    if (relation == GVA_EQUIVALENT || relation == GVA_DISJOINT)
                    {
                        relation = GVA_EQUIVALENT;
                    } // if
                } // if
                else if (node_parts_table[hash_idx].relation == GVA_CONTAINS)
                {
                    included += node_parts_table[hash_idx].included;
                    if (relation == GVA_IS_CONTAINED)
                    {
                        relation = GVA_OVERLAP;  // TODO: included = 1; ???
                        break;
                    } // if
                    relation = GVA_CONTAINS;
                } // if
                else if (node_parts_table[hash_idx].relation == GVA_IS_CONTAINED)
                {
                    size_t const part_idx = node_parts_table[hash_idx].start;
                    if (relation == GVA_CONTAINS)
                    {
                        included += node_parts_table[hash_idx].included;  // TODO: included = 1; ???
                        relation = GVA_OVERLAP;
                        break;
                    } // if

                    if (part_idx != is_contained_part_idx)
                    {
                        repair_is_contained(gva_std_allocator, reference, rhs_graph, tree, trie,
                                            &is_contained_nodes, is_contained_part_idx, &included, &relation);
                        is_contained_part_idx = part_idx;
                    } // if
                    ARRAY_APPEND(gva_std_allocator, is_contained_nodes, node_idx);
                    relation = GVA_IS_CONTAINED;
                } // if
                else if (node_parts_table[hash_idx].relation == GVA_OVERLAP)
                {
                    included = 1;
                    relation = GVA_OVERLAP;
                    break;
                } // if
            } // for all nodes for this allele
            repair_is_contained(gva_std_allocator, reference, rhs_graph, tree, trie,
                                &is_contained_nodes, is_contained_part_idx, &included, &relation);

            if (included > 0)
            {
                gva_uint const lhs_excluded = db_alleles[allele_idx].distance - included;
                gva_uint const rhs_excluded = rhs_graph.distance - included;

                if (lhs_excluded > 0 && rhs_excluded > 0)
                {
                    relation = GVA_OVERLAP;
                } // if
                else if (lhs_excluded > 0)
                {
                    relation = GVA_CONTAINS;
                } // if
                else if (rhs_excluded > 0)
                {
                    relation = GVA_IS_CONTAINED;
                } // if
                else // if (lhs_excluded == 0 && rhs_excluded == 0)
                {
                    relation = GVA_EQUIVALENT;
                } // else

                // only for testing
                if (relation != GVA_EQUIVALENT || db_alleles[allele_idx].line < line_count)
                {
                    printf("%u %zu %s\n", db_alleles[allele_idx].line, line_count, GVA_RELATION_LABELS[relation]);
                } // if
            } // if
        } // for all alleles

        results_table = HASH_TABLE_DESTROY(gva_std_allocator, results_table);
        node_parts_table = HASH_TABLE_DESTROY(gva_std_allocator, node_parts_table);

        gva_string_destroy(gva_std_allocator, rhs_graph.observed);
        gva_lcs_graph_destroy(gva_std_allocator, rhs_graph);

        line_count += 1;
    } // while fgets
    fclose(stream);

    db_alleles = ARRAY_DESTROY(gva_std_allocator, db_alleles);
    node_allele_join = ARRAY_DESTROY(gva_std_allocator, node_allele_join);
    interval_tree_destroy(gva_std_allocator, &tree);
    trie_destroy(gva_std_allocator, &trie);
    gva_string_destroy(gva_std_allocator, reference);

    return EXIT_SUCCESS;
} // dbsnp_main


int
locals_main(int argc, char* argv[static argc + 1])
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
        gva_uint data;
        gva_uint join_start;  // offset into node_allele_join
        gva_uint distance;
    }* db_alleles = NULL;

    errno = 0;
    stream = fopen(argv[2], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    static char line[LINE_SIZE] = {0};
    size_t line_count = 0;
    size_t parts = 0;
    while (fgets(line, sizeof(line), stream) != NULL)
    {
        size_t idx = 0;
        size_t data = parse_number(line, &idx);  // glob index
        idx += 1;  // skip space or tab
        size_t allele_entries = parse_number(line, &idx);
        parts += allele_entries;

        size_t const allele_idx = ARRAY_APPEND(gva_std_allocator, db_alleles, ((struct Allele) {data, array_length(node_allele_join), 0})) - 1;
        fprintf(stderr, "allele: %zu (%zu) length: %zu\n", data, allele_idx, allele_entries);

        size_t allele_distance = 0;
        line_count += 1;
        for (size_t i = 0; i < allele_entries; ++i)
        {
            fgets(line, sizeof(line), stream);
            idx = 0;
            size_t const distance = parse_number(line, &idx);
            allele_distance += distance;
            idx += 1;
            int const len = (char*) memchr(line + idx, '\n', LINE_SIZE - idx) - (line + idx);
            GVA_Variant variant;
            if (gva_parse_spdi(len, line + idx, &variant) == 0) {
                fprintf(stderr, "error: SPDI parsing failed at line %zu: %s", line_count + 1, line);
                continue;
            } //
            // fprintf(stderr, GVA_VARIANT_FMT " %zu\n", GVA_VARIANT_PRINT(variant), distance);

            gva_uint const inserted_idx = trie_insert(gva_std_allocator, &trie, variant.sequence.len, variant.sequence.str);
            gva_uint const tmp_idx = ARRAY_APPEND(gva_std_allocator, tree.nodes,
                                                  ((Interval_Tree_Node) {{GVA_NULL, GVA_NULL},
                                                  variant.start, variant.end, variant.end, 0, inserted_idx, GVA_NULL, distance})) - 1;
            gva_uint const node_idx = interval_tree_insert(&tree, tmp_idx);
            if (node_idx != tmp_idx)
            {
                array_header(tree.nodes)->length -= 1;  // reverse; node already in the tree
            } // if
            tree.nodes[node_idx].alleles = ARRAY_APPEND(gva_std_allocator, node_allele_join, ((struct Node_Allele) {node_idx, allele_idx, tree.nodes[node_idx].alleles})) - 1;

            line_count += 1;
        } // for
    } // while
    fprintf(stderr, "line count: %zu\n", line_count);

    fprintf(stderr, "parts: %zu\n", parts);
    fprintf(stderr, "tree nodes: %zu\n", array_length(tree.nodes));
    fprintf(stderr, "trie nodes: %zu (%zu)\n", array_length(trie.nodes), trie.strings.len);

    fprintf(stderr, "#db_alleles: %zu\n", array_length(db_alleles));
    fprintf(stderr, "#join:    %zu\n", array_length(node_allele_join));

    // line_count = 0;
    // parts = 0;
    // while (fgets(line, sizeof(line), stdin) != NULL) {
    //     size_t idx = 0;
    //     size_t data = parse_number(line, &idx);  // glob index
    //     idx += 1;  // skip space or tab
    //     size_t allele_entries = parse_number(line, &idx);
    //     parts += allele_entries;
    //
    //     size_t const allele_idx = ARRAY_APPEND(gva_std_allocator, db_alleles,
    //                                            ((struct Allele) {data, array_length(node_allele_join), 0})) - 1;
    //     fprintf(stderr, "allele: %zu (%zu) length: %zu\n", data, allele_idx, allele_entries);
    //
    //     // Join nodes in the index to parts in the query
    //     struct NODE_PARTS
    //     {
    //         HASH_TABLE_KEY;
    //         GVA_Relation relation;
    //         gva_uint start;
    //         gva_uint end;
    //         gva_uint included;
    //     }* node_parts_table = hash_table_init(gva_std_allocator, 1024, sizeof(*node_parts_table));
    //
    //     gva_uint part_distances[allele_entries];
    //     gva_uint part_starts[allele_entries];
    //     gva_uint part_ends[allele_entries];
    //
    //     size_t allele_distance = 0;
    //     line_count += 1;
    //     for (size_t part_idx = 0; part_idx < allele_entries; ++part_idx) {
    //         fgets(line, sizeof(line), stdin);
    //         idx = 0;
    //         size_t const rhs_distance = parse_number(line, &idx);
    //         allele_distance += rhs_distance;
    //         idx += 1;
    //         int const len = (char *) memchr(line + idx, '\n', LINE_SIZE - idx) - (line + idx);
    //         GVA_Variant rhs_part;
    //         if (gva_parse_spdi(len, line + idx, &rhs_part) == 0) {
    //             fprintf(stderr, "error: SPDI parsing failed at line %zu: %s", line_count + 1, line);
    //             continue;
    //         } //
    //         // fprintf(stderr, GVA_VARIANT_FMT " %zu\n", GVA_VARIANT_PRINT(rhs_part), distance);
    //         part_distances[part_idx] = rhs_distance;
    //         part_starts[part_idx] = rhs_part.start;
    //         part_ends[part_idx] = rhs_part.end;
    //
    //         gva_uint* candidates = interval_tree_intersection(gva_std_allocator, tree, rhs_part.start, rhs_part.end);
    //         for (size_t can_idx = 0; can_idx < array_length(candidates); ++can_idx)
    //         {
    //             gva_uint const node_idx = candidates[can_idx];
    //             GVA_Variant const db_var = {tree.nodes[node_idx].start, tree.nodes[node_idx].end,
    //                                         trie_string(trie, tree.nodes[node_idx].inserted)};
    //             GVA_Relation const relation = compare_from_index(reference, db_var, tree.nodes[node_idx].distance, rhs_part, rhs_distance);
    //             if (relation == GVA_DISJOINT)
    //             {
    //                 continue;
    //             } // if
    //
    //             // link nodes to query parts
    //             size_t hash_idx = HASH_TABLE_INDEX(node_parts_table, node_idx);
    //             if (node_parts_table[hash_idx].gva_key != node_idx)
    //             {
    //                 HASH_TABLE_SET(gva_std_allocator, node_parts_table, node_idx,
    //                                ((struct NODE_PARTS) {node_idx, relation, part_idx, part_idx + 1, rhs_distance}));
    //                 hash_idx = HASH_TABLE_INDEX(node_parts_table, node_idx);
    //             } // if
    //
    //             if (relation == GVA_EQUIVALENT || relation == GVA_IS_CONTAINED)
    //             {
    //                 node_parts_table[hash_idx].included = tree.nodes[node_idx].distance;
    //             } // if
    //             else if (relation == GVA_CONTAINS)
    //             {
    //                 node_parts_table[hash_idx].end = part_idx + 1;
    //             } // if
    //             else // GVA_OVERLAP
    //             {
    //                 node_parts_table[hash_idx].included = 1;
    //             } // else
    //         } // for every candidate
    //         candidates = ARRAY_DESTROY(gva_std_allocator, candidates);
    //
    //     } // for
    //
    //     // fix containment for multiple parts in single node for every query
    //     for (size_t npt_index = 0; npt_index < array_header(node_parts_table)->capacity; ++npt_index)
    //     {
    //         size_t const node_idx = node_parts_table[npt_index].gva_key;
    //         if (node_idx == GVA_NOT_FOUND)
    //         {
    //             continue;
    //         } // if
    //
    //         if (node_parts_table[npt_index].relation == GVA_CONTAINS &&
    //             node_parts_table[npt_index].end - node_parts_table[npt_index].start > 1)
    //         {
    //             size_t const lhs_distance = tree.nodes[node_idx].distance;
    //             size_t rhs_distance = 0;
    //             for (size_t i = node_parts_table[npt_index].start; i < node_parts_table[npt_index].end; ++i)
    //             {
    //                 // TODO: which distance?!
    //                 // rhs_distance += rhs_graph.local_supremal[i + 1].distance;
    //                 rhs_distance += part_distances[i];
    //
    //             } // for
    //             if (rhs_distance >= lhs_distance)
    //             {
    //                 node_parts_table[npt_index].included = 1;
    //                 node_parts_table[npt_index].relation = GVA_OVERLAP;
    //                 continue;
    //             } // if
    //
    //             GVA_Variant const lhs = {tree.nodes[node_idx].start,
    //                                      tree.nodes[node_idx].end,
    //                                      trie_string(trie, tree.nodes[node_idx].inserted)};
    //
    //             GVA_Variant rhs;
    //             gva_edges(rhs_graph.observed.str,
    //                       // rhs_graph.local_supremal[node_parts_table[npt_index].start],
    //                       // rhs_graph.local_supremal[node_parts_table[npt_index].end],
    //                       rhs_graph.local_supremal[node_parts_table[npt_index].start],
    //                       rhs_graph.local_supremal[node_parts_table[npt_index].end],
    //                       node_parts_table[npt_index].start == 0,
    //                       node_parts_table[npt_index].end == array_length(rhs_graph.local_supremal) - 1,
    //                       &rhs);
    //
    //             size_t const distance = variants_distance(gva_std_allocator, reference.len, reference.str, lhs, rhs);
    //             if (lhs_distance - distance == rhs_distance)
    //             {
    //                 node_parts_table[npt_index].included = rhs_distance;
    //             } // if
    //             else
    //             {
    //                 node_parts_table[npt_index].included = 1;
    //                 node_parts_table[npt_index].relation = GVA_OVERLAP;
    //             } // if
    //         } // if
    //     } // for node_parts_table
    //
    //     node_parts_table = HASH_TABLE_DESTROY(gva_std_allocator, node_parts_table);
    // } // while




    db_alleles = ARRAY_DESTROY(gva_std_allocator, db_alleles);
    node_allele_join = ARRAY_DESTROY(gva_std_allocator, node_allele_join);
    interval_tree_destroy(gva_std_allocator, &tree);
    trie_destroy(gva_std_allocator, &trie);
    gva_string_destroy(gva_std_allocator, reference);

    return EXIT_SUCCESS;
} // dbsnp_main


int
main(int argc, char* argv[static argc + 1])
{
    // return fasta_blob_read(argc, argv);
    // return fasta_blob_write(argc, argv);
    // return vcf_main(argc, argv);
    return dbsnp_main(argc, argv);
    // return locals_main(argc, argv);
} // main
