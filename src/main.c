#include <errno.h>      // errno
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint8_t
#include <stdio.h>      // FILE, stderr, stdout, fclose, fopen, fprintf
#include <stdlib.h>     // EXIT_*, atoll, qsort
#include <string.h>     // strerror, strlen

#include "../include/compare.h"     // gva_compare_distance
#include "../include/index.h"       // GVA_Index, gva_index_*, GVA_Query_Result
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_lcs_graph_*
#include "../include/relations.h"   // GVA_RELATION_LABELS
#include "../include/parser.h"      // gva_parse_spdi
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/string.h"      // GVA_String, gva_string_destroy
#include "../include/utils.h"       // gva_fasta_sequence*
#include "../include/variant.h"     // GVA_VARIANT_*, GVA_Variant, gva_patch, gva_variant_*
#include "array.h"          // ARRAY_*, array_length
#include "common.h"         // MAX, MIN
#include "dfa.h"            // dfa_*
#include "trie.h"           // Trie, trie_*


#include <assert.h>     // DEBUG


#define LINE_SIZE 4096

// #define REFERENCE_ID "NC_000022.11"
// #define REFERENCE_ID "NC_000006.12"
#define REFERENCE_ID "NC_000001.11"


static inline size_t
repeats(GVA_Allocator const allocator,
    size_t const len, char const word[static len])
{
    size_t* lps = allocator.allocate(allocator.context, NULL, 0, len * sizeof(*lps));
    if (lps == NULL)
    {
        return 0;  // OOM
    } // if

    lps[0] = 0;
    size_t length = 0;
    size_t idx = 1;
    while (idx < len)
    {
        if (word[idx] == word[length])
        {
            length += 1;
            lps[idx] = length;
            idx += 1;
        } // if
        else if (length > 0)
        {
            length = lps[length - 1];
        } // if
        else
        {
            lps[idx] = 0;
            idx += 1;
        } // else
    } // while

    lps = allocator.allocate(allocator.context, lps, len * sizeof(*lps), 0);
    return len - length;
} // repeats


static inline size_t
hgvs_location(size_t const len, char buffer[static len],
    size_t const start, size_t const end)
{
    if (end - start == 0)
    {
        return snprintf(buffer, len, "%zu_%zu", start, start + 1);
    } // if
    if (end - start == 1)
    {
        return snprintf(buffer, len, "%zu", start + 1);
    } // if
    return snprintf(buffer, len, "%zu_%zu", start + 1, end);
} // hgvs_location


static size_t
hgvs_variant(size_t const len, char buffer[static len],
    size_t const len_ref, char const reference[static restrict len_ref],
    GVA_Variant const variant)
{
    size_t const inserted = repeats(gva_std_allocator, variant.sequence.len, variant.sequence.str);
    size_t const deleted = repeats(gva_std_allocator, variant.end - variant.start, reference + variant.start);

    fprintf(stderr, "%zu vs %zu\n", inserted, deleted);


    size_t idx = hgvs_location(len, buffer, variant.start, variant.end);

    if (variant.end - variant.start == 0)
    {
        if (variant.sequence.len == 0)
        {
            return idx + snprintf(buffer + idx, len - idx, "=");
        } // if
        return idx + snprintf(buffer + idx, len - idx, "ins" GVA_STRING_FMT, GVA_STRING_PRINT(variant.sequence));
    } // if

    if (variant.end - variant.start == 1)
    {
        if (variant.sequence.len == 0)
        {
            return idx + snprintf(buffer + idx, len - idx, "del");
        } // if
        if (variant.sequence.len == 1)
        {
            return idx + snprintf(buffer + idx, len - idx, "%c>" GVA_STRING_FMT, reference[variant.start], GVA_STRING_PRINT(variant.sequence));
        } // if
        return idx + snprintf(buffer + idx, len - idx, "delins" GVA_STRING_FMT, GVA_STRING_PRINT(variant.sequence));
    } // if

    if (variant.sequence.len == 0)
    {
        return idx + snprintf(buffer + idx, len - idx, "del");
    } // if

    return idx + snprintf(buffer + idx, len - idx, "delins" GVA_STRING_FMT, GVA_STRING_PRINT(variant.sequence));
} // hgvs_variant


static size_t
to_hgvs(size_t const len, char buffer[static len],
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const n, GVA_Variant const variants[static restrict n])
{
    if (n == 0)
    {
        return snprintf(buffer, len, "=");
    } // if

    if (n == 1)
    {
        return hgvs_variant(len, buffer, len_ref, reference, variants[0]);
    } // if

    size_t idx = snprintf(buffer, len, "[");
    for (size_t i = 0; i < n; ++i)
    {
        idx += hgvs_variant(len - idx, buffer + idx, len_ref, reference, variants[i]);
        if (i < n - 1)
        {
            idx += snprintf(buffer + idx, len - idx, ";");
        } // if
    } // for
    return idx + snprintf(buffer + idx, len - idx, "]");
} // to_hgvs


// line: alphanumeric_id SPDI [distance]
static bool
parse_line(char const line[static restrict LINE_SIZE],
    GVA_String id[static restrict 1], GVA_Variant variant[static restrict 1], size_t distance[static restrict 1])
{
    id->len = strcspn(line, "\t ");
    if (id->len == 0)
    {
        return false;
    } // if
    id->str = line;

    size_t const len_spdi = strcspn(line + id->len + 1, "\n\t ");
    if (gva_parse_spdi(len_spdi, line + id->len + 1, variant) == 0)
    {
        return false;
    } // if

    *distance = atoll(line + id->len + len_spdi + 2);
    return true;
} // parse_line


static int
compare_query_alleles(void const* a, void const* b)
{
    GVA_Query_Allele const lhs = *(GVA_Query_Allele*) a;
    GVA_Query_Allele const rhs = *(GVA_Query_Allele*) b;

    if (lhs.included > rhs.included ||
        (lhs.included == rhs.included && lhs.excluded < rhs.excluded))
    {
        return -1;
    } // if
    return 1;
} // compare_query_alleles


static void
report_query_result(GVA_Index const* const index, GVA_LCS_Graph const graph, GVA_String const id,
    GVA_Query_Result const result,
    bool const hits, bool const query_disjoint, bool const index_disjoint)
{
    gva_uint* join = NULL;
    gva_uint* query = NULL;

    for (size_t i = 0; i < array_length(result.alleles); ++i)
    {
        fprintf(stdout, GVA_STRING_FMT " %s " GVA_STRING_FMT " %u %u %zu\n",
            GVA_STRING_PRINT(gva_index_allele_id(index, result.alleles[i].idx)),
            GVA_RELATION_LABELS[result.alleles[i].relation],
            GVA_STRING_PRINT(id),
            result.alleles[i].included, result.alleles[i].excluded, array_length(graph.dom_nodes) - 1);

        if (hits)
        {
            fprintf(stderr, GVA_STRING_FMT " " GVA_STRING_FMT " %u %u %s\n",
                GVA_STRING_PRINT(gva_index_allele_id(index, result.alleles[i].idx)),
                GVA_STRING_PRINT(id),
                result.alleles[i].included, result.alleles[i].excluded,
                GVA_RELATION_LABELS[result.alleles[i].relation]);

            gva_uint prev_join = gva_index_allele_parts(index, result.alleles[i].idx).start;
            size_t distance_join = 0;
            gva_uint prev_query = 0;
            size_t distance_query = 0;
            if (join != NULL)
            {
                array_header(join)->length = 0;
            } // if
            if (query != NULL)
            {
                array_header(query)->length = 0;
            } // if

            for (size_t j = result.alleles[i].hits.start; j < result.alleles[i].hits.end; ++j)
            {
                if (result.hits[j].index_parts.end - result.hits[j].index_parts.start > 1)
                {
                    fprintf(stderr, "  " REFERENCE_ID ":[");
                    size_t sum = 0;
                    for (gva_uint k = result.hits[j].index_parts.start; k < result.hits[j].index_parts.end; ++k)
                    {
                        if (index_disjoint)
                        {
                            for (gva_uint ii = prev_join; ii < result.parts[k]; ++ii)
                            {
                                ARRAY_APPEND(gva_std_allocator, join, ii);
                                distance_join += gva_index_variant_distance(index, result.alleles[i].idx, ii);
                            } // for
                            prev_join = result.parts[k] + 1;
                        } // if
                        fprintf(stderr, GVA_VARIANT_FMT_SPDI_ALLELE " (%zu), ",
                            GVA_VARIANT_PRINT_SPDI_ALLELE(gva_index_variant(index, result.alleles[i].idx, result.parts[k])),
                            gva_index_variant_distance(index, result.alleles[i].idx, result.parts[k]));
                        sum += gva_index_variant_distance(index, result.alleles[i].idx, result.parts[k]);
                    } // for
                    fprintf(stderr, "] (%zu) " GVA_VARIANT_FMT_SPDI " (%u) %u %u %s\n",
                        sum,
                        GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, gva_lcs_graph_local_supremal(graph, result.hits[j].query_parts.start, result.hits[j].query_parts.start + 1)),
                        graph.dom_nodes[result.hits[j].query_parts.start + 1].distance - graph.dom_nodes[result.hits[j].query_parts.start].distance,
                        result.hits[j].included, result.hits[j].excluded, GVA_RELATION_LABELS[result.hits[j].relation]);
                    if (query_disjoint)
                    {
                        for (gva_uint k = prev_query; k < result.hits[j].query_parts.start; ++k)
                        {
                            ARRAY_APPEND(gva_std_allocator, query, k);
                        } // for
                        distance_query += graph.dom_nodes[result.hits[j].query_parts.start].distance - graph.dom_nodes[prev_query].distance;
                        prev_query = result.hits[j].query_parts.start + 1;
                    } // if
                } // if
                else if (result.hits[j].query_parts.end - result.hits[j].query_parts.start > 1)
                {
                    fprintf(stderr, "  " GVA_VARIANT_FMT_SPDI " (%zu) " REFERENCE_ID ":[",
                        GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, gva_index_variant(index, result.alleles[i].idx, result.hits[j].index_parts.start)),
                        gva_index_variant_distance(index, result.alleles[i].idx, result.hits[j].index_parts.start));
                    size_t sum = 0;
                    for (gva_uint k = result.hits[j].query_parts.start; k < result.hits[j].query_parts.end; ++k)
                    {
                        if (query_disjoint)
                        {
                            for (gva_uint ii = prev_query; ii < result.parts[k]; ++ii)
                            {
                                ARRAY_APPEND(gva_std_allocator, query, ii);
                            } // for
                            distance_query += graph.dom_nodes[result.parts[k]].distance - graph.dom_nodes[prev_query].distance;
                            prev_query = result.parts[k] + 1;
                        } // if
                        fprintf(stderr, GVA_VARIANT_FMT_SPDI_ALLELE " (%u), ",
                            GVA_VARIANT_PRINT_SPDI_ALLELE(gva_lcs_graph_local_supremal(graph, result.parts[k], result.parts[k] + 1)),
                            graph.dom_nodes[result.parts[k] + 1].distance - graph.dom_nodes[result.parts[k]].distance);
                        sum += graph.dom_nodes[result.parts[k] + 1].distance - graph.dom_nodes[result.parts[k]].distance;
                    } // for
                    fprintf(stderr, "] (%zu) %u %u %s\n", sum, result.hits[j].included, result.hits[j].excluded, GVA_RELATION_LABELS[result.hits[j].relation]);
                    if (index_disjoint)
                    {
                        for (gva_uint k = prev_join; k < result.hits[j].index_parts.start; ++k)
                        {
                            ARRAY_APPEND(gva_std_allocator, join, k);
                            distance_join += gva_index_variant_distance(index, result.alleles[i].idx, k);
                        } // for
                        prev_join = result.hits[j].index_parts.start + 1;
                    } // if
                } // if
                else if (result.hits[j].included > 0)
                {
                    fprintf(stderr, "  " GVA_VARIANT_FMT_SPDI " (%u) " GVA_VARIANT_FMT_SPDI " (%zu) %u %u %s\n",
                        GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, gva_lcs_graph_local_supremal(graph, result.hits[j].query_parts.start, result.hits[j].query_parts.start + 1)),
                        graph.dom_nodes[result.hits[j].query_parts.start + 1].distance - graph.dom_nodes[result.hits[j].query_parts.start].distance,
                        GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, gva_index_variant(index, result.alleles[i].idx, result.hits[j].index_parts.start)),
                        gva_index_variant_distance(index, result.alleles[i].idx, result.hits[j].index_parts.start),
                        result.hits[j].included, result.hits[j].excluded, GVA_RELATION_LABELS[result.hits[j].relation]);
                    if (index_disjoint)
                    {
                        for (gva_uint k = prev_join; k < result.hits[j].index_parts.start; ++k)
                        {
                            ARRAY_APPEND(gva_std_allocator, join, k);
                            distance_join += gva_index_variant_distance(index, result.alleles[i].idx, k);
                        } // for
                        prev_join = result.hits[j].index_parts.start + 1;
                    } // if
                    if (query_disjoint)
                    {
                        for (gva_uint k = prev_query; k < result.hits[j].query_parts.start; ++k)
                        {
                            ARRAY_APPEND(gva_std_allocator, query, k);
                        } // for
                        distance_query += graph.dom_nodes[result.hits[j].query_parts.start].distance - graph.dom_nodes[prev_query].distance;
                        prev_query = result.hits[j].query_parts.start + 1;
                    } // if
                } // if
            } // for

            if (index_disjoint)
            {
                for (gva_uint k = prev_join; k < gva_index_allele_parts(index, result.alleles[i].idx).end; ++k)
                {
                    ARRAY_APPEND(gva_std_allocator, join, k);
                    distance_join += gva_index_variant_distance(index, result.alleles[i].idx, k);
                } // for
            } // if
            if (query_disjoint)
            {
                for (gva_uint k = prev_query; k < array_length(graph.dom_nodes) - 1; ++k)
                {
                    ARRAY_APPEND(gva_std_allocator, query, k);
                } // for
                distance_query += graph.dom_nodes[array_length(graph.dom_nodes) - 1].distance - graph.dom_nodes[prev_query].distance;
            } // if

            if (index_disjoint && join != NULL)
            {
                fprintf(stderr, "  " REFERENCE_ID ":[");
                for (size_t j = 0; j < array_length(join); ++j)
                {
                    fprintf(stderr, GVA_VARIANT_FMT_SPDI_ALLELE " (%zu), ",
                        GVA_VARIANT_PRINT_SPDI_ALLELE(gva_index_variant(index, result.alleles[i].idx, join[j])),
                        gva_index_variant_distance(index, result.alleles[i].idx, join[j]));
                } // for
                fprintf(stderr, "] * 0 %zu disjoint\n", distance_join);
            } // if
            if (query_disjoint && query != NULL)
            {
                fprintf(stderr, "  * " REFERENCE_ID ":[");
                for (size_t j = 0; j < array_length(query); ++j)
                {
                    fprintf(stderr, GVA_VARIANT_FMT_SPDI_ALLELE " (%u), ",
                        GVA_VARIANT_PRINT_SPDI_ALLELE(gva_lcs_graph_local_supremal(graph, query[j], query[j] + 1)),
                        graph.dom_nodes[query[j] + 1].distance - graph.dom_nodes[query[j]].distance);
                } // for
                fprintf(stderr, "] 0 %zu disjoint\n", distance_query);
            } // if
        } // if
    } // for

    join = ARRAY_DESTROY(gva_std_allocator, join);
    query = ARRAY_DESTROY(gva_std_allocator, query);
} // report_query_result


int
index_main(int argc, char* argv[static argc])
{
    if (argc < 3)
    {
        fprintf(stderr, "usage %s reference.blob data\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    errno = 0;
    FILE* stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    GVA_String reference = {0};
    reference = gva_fasta_sequence_blob(gva_std_allocator, stream);
    fclose(stream);

    fprintf(stderr, "reference length: %zu\n", reference.len);

    GVA_Index* index = gva_index_init(gva_std_allocator, reference.len, reference.str);
    if (index == NULL)
    {
        fprintf(stderr, "OOM\n");
        gva_string_destroy(gva_std_allocator, reference);
        return EXIT_FAILURE;
    } // if

    errno = 0;
    stream = fopen(argv[2], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        gva_string_destroy(gva_std_allocator, reference);
        return EXIT_FAILURE;
    } // if

    size_t line_count = 0;
    static char line[LINE_SIZE] = {0};
    while (fgets(line, sizeof(line), stream) != NULL)
    {
        line_count += 1;
        GVA_String id = {0};
        GVA_Variant variant = {0};
        size_t distance = 0;
        if (!parse_line(line, &id, &variant, &distance))
        {
            fprintf(stderr, "parsing failed at line %zu: %s", line_count, line);
            continue;
        } // if
        gva_index_insert(index, id.len, id.str, variant, distance);
    } // while
    fclose(stream);

    fprintf(stderr, "Index populated\n");
    gva_index_stats(index);

    line_count = 0;
    while (fgets(line, sizeof(line), stdin) != NULL)
    {
        line_count += 1;
        GVA_String id = {0};
        GVA_Variant variant = {0};
        if (!parse_line(line, &id, &variant, &(size_t) {0}))
        {
            fprintf(stderr, "parsing failed at line %zu: %s", line_count, line);
            continue;
        } // if
        GVA_LCS_Graph graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, 1, &variant);

        //fprintf(stderr, "\nQuery (" GVA_STRING_FMT "): " GVA_VARIANT_FMT_SPDI " (%zu)\n", GVA_STRING_PRINT(id), GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, gva_lcs_graph_supremal(graph)), gva_lcs_graph_distance(graph));
        GVA_Query_Result result = gva_index_query(gva_std_allocator, index, graph);

        if (array_length(result.alleles) > 0)
        {
            qsort(result.alleles, array_length(result.alleles), sizeof(*result.alleles), compare_query_alleles);
        } // if

        report_query_result(index, graph, id, result, false, true, true);

        result.alleles = ARRAY_DESTROY(gva_std_allocator, result.alleles);
        result.hits = ARRAY_DESTROY(gva_std_allocator, result.hits);
        result.parts = ARRAY_DESTROY(gva_std_allocator, result.parts);
        gva_lcs_graph_destroy(gva_std_allocator, graph, true);
    } // while

    index = gva_index_destroy(index);
    gva_string_destroy(gva_std_allocator, reference);

    return EXIT_SUCCESS;
} // index_main


int
allele_main(int argc, char* argv[static argc])
{
    if (argc < 3)
    {
        fprintf(stderr, "usage: %s reference observed\n", argv[0]);
        return EXIT_FAILURE;
    } // if

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

    errno = 0;
    stream = fopen(argv[2], "r");
    if (stream == NULL)
    {
        gva_string_destroy(gva_std_allocator, reference);
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    GVA_Variant* variants = NULL;

    size_t line_count = 0;
    size_t dropped = 0;
    static char line[LINE_SIZE] = {0};
    while (fgets(line, sizeof(line), stream) != NULL)
    {
        line_count += 1;

        GVA_Variant variant = {0};
        if (gva_parse_spdi(strlen(line) - 1, line, &variant) == 0)
        {
            fprintf(stderr, "error: SPDI parsing failed at line %zu: %s", line_count, line);
            continue;
        } // if

        GVA_Variant const trimmed = gva_variant_prefix_trimmed(reference.len, reference.str, variant);

        if (array_length(variants) > 0 && trimmed.start < variants[array_length(variants) - 1].end)
        {
            fprintf(stderr, "dropped at line %zu: %s", line_count + 1, line);
            dropped += 1;
            continue;
        } // if

        ARRAY_APPEND(gva_std_allocator, variants, gva_variant_dup(gva_std_allocator, trimmed));
    } // while
    fclose(stream);
    fprintf(stderr, "#variants: %zu\n", array_length(variants));
    fprintf(stderr, "#dropped:  %zu\n", dropped);

    GVA_LCS_Graph graph = gva_lcs_graph_from_allele(gva_std_allocator, reference.len, reference.str, array_length(variants), variants);
    for (size_t i = 0; i < array_length(graph.dom_nodes) - 1; ++i)
    {
        fprintf(stdout, GVA_VARIANT_FMT_SPDI " %u\n",
            GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, gva_lcs_graph_local_supremal(graph, i, i + 1)),
            graph.dom_nodes[i + 1].distance - graph.dom_nodes[i].distance);
    } // for
    fprintf(stderr, "#nodes: %zu\n", array_length(graph.nodes));
    fprintf(stderr, "#edges: %zu\n", array_length(graph.edges));

    gva_lcs_graph_destroy(gva_std_allocator, graph, true);

    if (variants != NULL)
    {
        for (size_t i = 0; i < array_length(variants); ++i)
        {
            gva_string_destroy(gva_std_allocator, variants[i].sequence);
        } // for
        variants = ARRAY_DESTROY(gva_std_allocator, variants);
    } // if

    gva_string_destroy(gva_std_allocator, reference);

    return EXIT_SUCCESS;
} // allele_main


typedef struct
{
    gva_uint start;
    gva_uint end;
    gva_uint inserted;
    gva_uint distance;
    gva_uint label;
    gva_uint dfa;
} Entry;


static inline int
compare_entries(void const* a, void const* b)
{
    Entry const* const lhs = a;
    Entry const* const rhs = b;

    if (lhs->start < rhs->start)
    {
        return -1;
    } // if
    if (lhs->start > rhs->start)
    {
        return 1;
    } // if
    if (lhs->distance > rhs->distance)
    {
        return -1;
    } // if
    if (lhs->distance < rhs->distance)
    {
        return 1;
    } // if
    return 0;
} // compare_entries


int
all_main(int argc, char* argv[static argc])
{
    if (argc < 3)
    {
        fprintf(stderr, "usage %s reference.blob data\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    errno = 0;
    FILE* stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    GVA_String reference = {0};
    reference = gva_fasta_sequence_blob(gva_std_allocator, stream);
    fclose(stream);

    fprintf(stderr, "reference length: %zu\n", reference.len);

    errno = 0;
    stream = fopen(argv[2], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        gva_string_destroy(gva_std_allocator, reference);
        return EXIT_FAILURE;
    } // if

    Trie labels = trie_init(gva_std_allocator, gva_std_allocator, GVA_NULL);
    Trie sequences = trie_init(gva_std_allocator, gva_std_allocator, GVA_NULL);
    Trie dfas = trie_init(gva_std_allocator, gva_std_allocator, GVA_NULL);
    Entry* entries = NULL;

    size_t line_count = 0;
    static char line[LINE_SIZE] = {0};
    while (fgets(line, sizeof(line), stream) != NULL)
    {
        line_count += 1;
        GVA_String id = {0};
        GVA_Variant variant = {0};
        size_t distance = 0;
        if (!parse_line(line, &id, &variant, &distance))
        {
            fprintf(stderr, "parsing failed at line %zu: %s", line_count, line);
            continue;
        } // if

        if (distance == 0)
        {
            fprintf(stderr, "ERROR: DISTANCE\n");
            continue;
        } // if

        uint8_t* dfa = dfa_from_alignment(gva_std_allocator, NULL, variant.end - variant.start, reference.str + variant.start, variant.sequence.len, variant.sequence.str);

        ARRAY_APPEND(gva_std_allocator, entries,
            ((Entry)
            {
                .start = variant.start,
                .end = variant.end,
                .inserted = trie_insert(&sequences, variant.sequence.len, variant.sequence.str),
                .distance = distance,
                .label = trie_insert(&labels, id.len, id.str),
                .dfa = trie_insert(&dfas, array_length(dfa), (char*) dfa),
            }));

        dfa = ARRAY_DESTROY(gva_std_allocator, dfa);
    } // while
    fclose(stream);

    fprintf(stderr, "#lines: %zu\n", line_count);

    if (entries != NULL)
    {
        qsort(entries, array_length(entries), sizeof(*entries), compare_entries);
        fprintf(stderr, "#entries: %zu\n", array_length(entries));

        size_t count = 0;
        for (size_t i = 0; i < array_length(entries); ++i)
        {
            GVA_Variant const lhs = { entries[i].start, entries[i].end, trie_string(sequences, entries[i].inserted) };

            for (size_t j = i + 1; j < array_length(entries); ++j)
            {
                if (entries[i].end < entries[j].start)
                {
                    break;
                } // if
                count += 1;

                if (entries[i].start == entries[j].start && entries[i].end == entries[j].end && entries[i].inserted == entries[j].inserted)
                {
                    fprintf(stdout, GVA_STRING_FMT " equivalent " GVA_STRING_FMT "\n",
                        GVA_STRING_PRINT(trie_string(labels, entries[i].label)),
                        GVA_STRING_PRINT(trie_string(labels, entries[j].label)));
                    continue;  // equivalent
                } // if

                GVA_Variant const rhs = { entries[j].start, entries[j].end, trie_string(sequences, entries[j].inserted) };

                size_t const distance = gva_compare_distance(gva_std_allocator, reference.len, reference.str, lhs, rhs);

                if (entries[i].distance + entries[j].distance == distance)
                {
                    continue;  // disjoint
                } // if

                if (entries[i].distance - entries[j].distance == distance)
                {
                    fprintf(stdout, GVA_STRING_FMT " contains " GVA_STRING_FMT "\n",
                        GVA_STRING_PRINT(trie_string(labels, entries[i].label)),
                        GVA_STRING_PRINT(trie_string(labels, entries[j].label)));
                    continue;  // contains
                } // if

                if (entries[j].distance - entries[i].distance == distance)
                {
                    fprintf(stdout, GVA_STRING_FMT " is_contained " GVA_STRING_FMT "\n",
                        GVA_STRING_PRINT(trie_string(labels, entries[i].label)),
                        GVA_STRING_PRINT(trie_string(labels, entries[j].label)));
                    continue;  // is_contained
                } // if

                GVA_String lhs_dfa = trie_string(dfas, entries[i].dfa);
                GVA_String rhs_dfa = trie_string(dfas, entries[j].dfa);

                size_t const overlap = dfa_overlap(gva_std_allocator, reference.len, reference.str, lhs, (uint8_t*) lhs_dfa.str, rhs, (uint8_t*) rhs_dfa.str);
                // bool const disjoint = dfa_disjoint(lhs, (uint8_t*) lhs_dfa.str, rhs, (uint8_t*) rhs_dfa.str);

                // if (disjoint)
                if (overlap == 0)
                {
                    fprintf(stdout, GVA_STRING_FMT " disjoint " GVA_STRING_FMT "\n",
                        GVA_STRING_PRINT(trie_string(labels, entries[i].label)),
                        GVA_STRING_PRINT(trie_string(labels, entries[j].label)));
                } // if
                else
                {
                    fprintf(stdout, GVA_STRING_FMT " overlap " GVA_STRING_FMT " %zu\n",
                        GVA_STRING_PRINT(trie_string(labels, entries[i].label)),
                        GVA_STRING_PRINT(trie_string(labels, entries[j].label)),
                        overlap);
                } // else
            } // for
        } // for

        fprintf(stderr, "#combos: %zu\n", count);
    } // if

    entries = ARRAY_DESTROY(gva_std_allocator, entries);
    trie_destroy(&dfas);
    trie_destroy(&labels);
    trie_destroy(&sequences);

    gva_string_destroy(gva_std_allocator, reference);

    return EXIT_SUCCESS;
} // all_main


int
hgvs_main(int argc, char* argv[static argc])
{
    static char buffer[1024] = {'\0'};
    GVA_String const ref = {10, "AAATAATATAATAATTTAT"};
    GVA_Variant const variant = {2, 13, {6, "ATAATA"}};
    fprintf(stderr, "%zu\n", to_hgvs(1024, buffer, ref.len, ref.str, 1, &variant));

    fprintf(stderr, "%s\n", buffer);

    return EXIT_SUCCESS;


    if (argc < 2)
    {
        fprintf(stderr, "usage %s reference.blob\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    errno = 0;
    FILE* stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    GVA_String reference = {0};
    reference = gva_fasta_sequence_blob(gva_std_allocator, stream);
    fclose(stream);

    fprintf(stderr, "reference length: %zu\n", reference.len);

    size_t line_count = 0;
    static char line[LINE_SIZE] = {0};
    while (fgets(line, sizeof(line), stdin) != NULL)
    {
        line_count += 1;

        size_t idx = 0;
        size_t tok = 0;

        tok = strcspn(line + idx, "\t ");
        idx += tok + 1;  // ignore identifier

        tok = strcspn(line + idx, "\t ");
        GVA_Variant variant = {0};
        if (gva_parse_spdi(tok, line + idx, &variant) == 0)
        {
            fprintf(stderr, "%zu: ERROR SPDI: %s", line_count, line);
            continue;
        } // if
        idx += tok + 1;

        tok = strcspn(line + idx, "\n\t ");
        GVA_HGVS_Allele allele = gva_parse_hgvs(gva_std_allocator, reference.len, reference.str, tok, line + idx );
        if (!allele.interpretable)
        {
            fprintf(stderr, "%zu: ERROR HGVS: %s", line_count, line);
            continue;
        } // if
        idx += tok + 1;

        ARRAY_DESTROY(gva_std_allocator, allele.inserted);
        ARRAY_DESTROY(gva_std_allocator, allele.variants);
    } //while

    gva_string_destroy(gva_std_allocator, reference);

    return EXIT_SUCCESS;
} // hgvs_main


int
main(int argc, char* argv[static argc])
{
    // return hgvs_main(argc, argv);
    // return allele_main(argc, argv);
    return index_main(argc, argv);
    // return all_main(argc, argv);
} // main
