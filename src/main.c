#include <errno.h>      // errno
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // FILE, stderr, stdout, fclose, fopen, fprintf
#include <stdlib.h>     // EXIT_*, atoll, qsort
#include <string.h>     // strerror, strlen

#include "../include/compare.h"     // gva_compare_graphs
#include "../include/edit.h"        // gva_edit_distance
#include "../include/index.h"       // GVA_Index, gva_index_*, GVA_Query_Result
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_lcs_graph_*
#include "../include/relations.h"   // GVA_RELATION_LABELS
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/string.h"      // GVA_String, gva_string_destroy
#include "../include/utils.h"       // gva_fasta_sequence*, gva_lcs_graph_dot
#include "../include/variant.h"     // GVA_VARIANT_*, GVA_Variant, gva_parse_spdi, gva_variant_*
#include "align.h"          // LCS_Matches, lcs_align_one
#include "array.h"          // ARRAY_*, array_length
#include "bitset.h"         // bitset_*
#include "common.h"         // MAX, MIN
#include "nfa.h"            // DFA, dfa_*
#include "trie.h"           // Trie, trie_*


#include <assert.h>


#define LINE_SIZE 8194

// #define REFERENCE_ID "NC_000022.11"
// #define REFERENCE_ID "NC_000006.12"
#define REFERENCE_ID "NC_000001.11"


static inline bool
variant_disjoint(size_t const lhs_count, GVA_Variant const lhs,
    size_t const rhs_count, GVA_Variant const rhs)
{
    if (rhs.start > lhs.end + lhs_count - 1 ||
        lhs.start > rhs.end + rhs_count - 1)
    {
        return true;
    } // if

    size_t lhs_offset = 0;
    size_t rhs_offset = 0;
    if (rhs.start > lhs.start)
    {
        lhs_offset = MIN(rhs.start - lhs.start, lhs_count - 1);
    } // if
    else if (lhs.start > rhs.start)
    {
        rhs_offset = MIN(lhs.start - rhs.start, rhs_count - 1);
    } // if

    size_t const start = MAX(lhs.start + lhs_offset, rhs.start + rhs_offset);
    size_t const end = MIN(lhs.end + lhs_offset, rhs.end + rhs_offset);

    if (end > start)
    {
        return false;
    } // if

    unsigned __int128 lhs_seq = 0;
    for (size_t i = 0; i < lhs.sequence.len; ++i)
    {
        lhs_seq |= 1 << (lhs.sequence.str[lhs_offset + i] & 0x7f);
    } // for
    unsigned __int128 rhs_seq = 0;
    for (size_t i = 0; i < rhs.sequence.len; ++i)
    {
        rhs_seq |= 1 << (rhs.sequence.str[rhs_offset + i] & 0x7f);
    } // for

    return (lhs_seq & rhs_seq) == 0;
} // variant_disjoint


static bool
is_disjoint(GVA_Allocator const allocator,
    GVA_LCS_Graph const lhs, GVA_LCS_Graph const rhs)
{
    size_t const rhs_nodes = array_length(rhs.nodes);
    size_t const length = array_length(lhs.nodes) * rhs_nodes;
    struct
    {
        gva_uint seen;
        gva_uint next;
    }* queue = allocator.allocate(allocator.context, NULL, 0, length * sizeof(*queue));
    if (queue == NULL)
    {
        return true;  // FIXME: OOM
    } // if

    for (size_t i = 0; i < length; ++i)
    {
        queue[i].seen = false;
        queue[i].next = GVA_NULL;
    } // for

    gva_uint tail = lhs.source * rhs_nodes + rhs.source;
    for (gva_uint head = tail; head != GVA_NULL; head = queue[head].next)
    {
        size_t const lhs_idx = head / rhs_nodes;
        size_t const rhs_idx = head % rhs_nodes;

        fprintf(stderr, "{%zu, %zu}\n", lhs_idx, rhs_idx);

        for (gva_uint i = lhs.nodes[lhs_idx].edges; i != GVA_NULL; i = lhs.edges[i].next)
        {
            GVA_Variant lhs_variant = {0};
            size_t const lhs_count = gva_edges(lhs.observed.str,
                lhs.nodes[lhs_idx].match, lhs.nodes[lhs.edges[i].tail].match,
                lhs_idx == lhs.source, lhs.nodes[lhs.edges[i].tail].edges == GVA_NULL,
                &lhs_variant);

            for (gva_uint j = rhs.nodes[rhs_idx].edges; j != GVA_NULL; j = rhs.edges[j].next)
            {
                GVA_Variant rhs_variant = {0};
                size_t const rhs_count = gva_edges(rhs.observed.str,
                    rhs.nodes[rhs_idx].match, rhs.nodes[rhs.edges[j].tail].match,
                    rhs_idx == rhs.source, rhs.nodes[rhs.edges[j].tail].edges == GVA_NULL,
                    &rhs_variant);

                //fprintf(stderr, "    " GVA_VARIANT_FMT " x %zu  vs  " GVA_VARIANT_FMT " x %zu\n",
                //    GVA_VARIANT_PRINT(lhs_variant), lhs_count, GVA_VARIANT_PRINT(rhs_variant), rhs_count);

                if (!variant_disjoint(lhs_count, lhs_variant, rhs_count, rhs_variant))
                {
                    queue = allocator.allocate(allocator.context, queue, length * sizeof(*queue), 0);
                    return false;
                } // if

                size_t const next = lhs.edges[i].tail * rhs_nodes + rhs.edges[j].tail;
                if (lhs.nodes[lhs.edges[i].tail].edges != GVA_NULL &&
                    rhs.nodes[rhs.edges[j].tail].edges != GVA_NULL &&
                    !queue[next].seen)
                {
                    queue[tail].next = next;
                    tail = next;
                    queue[next].seen = true;
                } // if
            } // for
        } // for

        if (lhs.nodes[lhs_idx].lambda != GVA_NULL)
        {
            size_t const next = lhs.nodes[lhs_idx].lambda * rhs_nodes + rhs_idx;
            if (!queue[next].seen)
            {
                //fprintf(stderr, "  lambda lhs: {%u, %zu}\n", lhs.nodes[lhs_idx].lambda, rhs_idx);
                queue[tail].next = next;
                tail = next;
                queue[next].seen = true;
            } // if
        } // if
        if (rhs.nodes[rhs_idx].lambda != GVA_NULL)
        {
            size_t const next = lhs_idx * rhs_nodes + rhs.nodes[rhs_idx].lambda;
            if (!queue[next].seen)
            {
                //fprintf(stderr, "  lambda rhs: {%zu, %u}\n", lhs_idx, rhs.nodes[rhs_idx].lambda);
                queue[tail].next = next;
                tail = next;
                queue[next].seen = true;
            } // if
        } // if

        for (gva_uint i = lhs.nodes[lhs_idx].edges; i != GVA_NULL; i = lhs.edges[i].next)
        {
            size_t const next = lhs.edges[i].tail * rhs_nodes + rhs_idx;
            if (lhs.nodes[lhs.edges[i].tail].edges != GVA_NULL && !queue[next].seen)
            {
                //fprintf(stderr, "  advance lhs: {%u, %zu}\n", lhs.edges[i].tail, rhs_idx);
                queue[tail].next = next;
                tail = next;
                queue[next].seen = true;
            } // if
        } // for
        for (gva_uint i = rhs.nodes[rhs_idx].edges; i != GVA_NULL; i = rhs.edges[i].next)
        {
            size_t const next = lhs_idx * rhs_nodes + rhs.edges[i].tail;
            if (rhs.nodes[rhs.edges[i].tail].edges != GVA_NULL && !queue[next].seen)
            {
                //fprintf(stderr, "  advance rhs: {%zu, %u}\n", lhs_idx, rhs.edges[i].tail);
                queue[tail].next = next;
                tail = next;
                queue[next].seen = true;
            } // if
        } // for

    } // for

    queue = allocator.allocate(allocator.context, queue, length * sizeof(*queue), 0);
    return true;
} // is_disjoint


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

        if (result.alleles[i].relation == GVA_OVERLAP &&
            result.alleles[i].hits.end - result.alleles[i].hits.start > 1)
        {
            bool overlap = false;
            for (size_t j = result.alleles[i].hits.start; j < result.alleles[i].hits.end; ++j)
            {
                if (result.hits[j].relation == GVA_OVERLAP)
                {
                    overlap = true;
                    break;
                } // if
            } // for
            if (!overlap)
            {
        fprintf(stdout, GVA_STRING_FMT " %s " GVA_STRING_FMT " %u %u %zu\n",
            GVA_STRING_PRINT(gva_index_allele_id(index, result.alleles[i].idx)),
            GVA_RELATION_LABELS[result.alleles[i].relation],
            GVA_STRING_PRINT(id),
            result.alleles[i].included, result.alleles[i].excluded, array_length(graph.dom_nodes) - 1);
            } // if
        } // if

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


int
supremal_main(int argc, char* argv[static argc])
{
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
        GVA_String id = {0};
        GVA_Variant variant = {0};
        size_t distance = 0;
        if (!parse_line(line, &id, &variant, &distance))
        {
            fprintf(stderr, "parsing failed at line %zu: %s", line_count, line);
            continue;
        } // if
        GVA_LCS_Graph graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, 1, &variant);

        GVA_Variant supremal = gva_lcs_graph_supremal(graph);

        fprintf(stdout, GVA_STRING_FMT " " GVA_VARIANT_FMT_SPDI " %zu\n",
            GVA_STRING_PRINT(id),
            GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, supremal),
            gva_lcs_graph_distance(graph));

        /*
        for (size_t i = 0; i < array_length(graph.dom_nodes) - 1; ++i)
        {
            fprintf(stdout, GVA_STRING_FMT " " GVA_VARIANT_FMT_SPDI " %u\n",
                GVA_STRING_PRINT(id),
                GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, gva_lcs_graph_local_supremal(graph, i, i + 1)),
                graph.dom_nodes[i + 1].distance - graph.dom_nodes[i].distance);
        } // for
        */

        gva_lcs_graph_destroy(gva_std_allocator, graph, true);
    } // while

    gva_string_destroy(gva_std_allocator, reference);

    return EXIT_SUCCESS;
} // supremal_main


int
overlap_main(int argc, char* argv[static argc])
{
    (void) argv;

    size_t line_count = 0;
    static char line[LINE_SIZE] = {0};
    while (fgets(line, sizeof(line), stdin) != NULL)
    {
        line_count += 1;
        size_t start = 0;
        size_t end = 0;

        end = strcspn(line, "\t ");
        GVA_String reference = {end, line + start};

        start += end + 1;
        end = strcspn(line + start, "\t ");

        GVA_Variant lhs_variant = {0};
        if (gva_parse_spdi(end, line + start, &lhs_variant) == 0)
        {
            fprintf(stderr, "Error parsing LHS SPDI at line %zu: %s", line_count, line);
            continue;
        } // if

        start += end + 1;
        end = strcspn(line + start, "\n\t ");

        GVA_Variant rhs_variant = {0};
        if (gva_parse_spdi(end, line + start, &rhs_variant) == 0)
        {
            fprintf(stderr, "Error parsing RHS SPDI at line %zu: %s", line_count, line);
            continue;
        } // if

        GVA_LCS_Graph lhs = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, 1, &lhs_variant);
        GVA_LCS_Graph rhs = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, 1, &rhs_variant);

        DFA lhs_dfa = dfa_from_lcs_graph(gva_std_allocator, lhs);
        DFA rhs_dfa = dfa_from_lcs_graph(gva_std_allocator, rhs);
        dfa_dot(lhs_dfa);
        dfa_dot(rhs_dfa);

        fprintf(stderr, "OVERLAP: %zu\n", dfa_max_overlap(gva_std_allocator, lhs_dfa, rhs_dfa));

        dfa_destroy(gva_std_allocator, &rhs_dfa);
        dfa_destroy(gva_std_allocator, &lhs_dfa);

        fprintf(stderr, "%s\n", GVA_RELATION_LABELS[gva_compare_graphs(gva_std_allocator, reference.len, reference.str, lhs, rhs)]);

        static char filename[128] = {0};
        snprintf(filename, 128, "overlap/graph_%zu.dot", line_count);
        FILE* stream = fopen(filename, "w");
        gva_lcs_graph_dot(stream, lhs);
        gva_lcs_graph_dot(stream, rhs);
        fclose(stream);

        fprintf(stderr, "%sdistances lhs: %zu rhs: %zu :: ", line, gva_lcs_graph_distance(lhs), gva_lcs_graph_distance(rhs));
        {
            size_t const distance = gva_edit_distance(gva_std_allocator,
                lhs_variant.sequence.len, lhs_variant.sequence.str,
                rhs_variant.sequence.len, rhs_variant.sequence.str);
            fprintf(stderr, "%zu\n", distance);
        }
        fprintf(stderr, "#nodes lhs: %zu rhs %zu :: %zu\n", array_length(lhs.nodes), array_length(rhs.nodes), array_length(lhs.nodes) * array_length(rhs.nodes));
        fprintf(stderr, "#edges lhs: %zu rhs %zu :: %zu\n", array_length(lhs.edges), array_length(rhs.edges), array_length(lhs.edges) * array_length(rhs.edges));

        //fprintf(stderr, "disjoint: %d\n", is_disjoint(gva_std_allocator, lhs, rhs));

        gva_lcs_graph_destroy(gva_std_allocator, rhs, true);
        gva_lcs_graph_destroy(gva_std_allocator, lhs, true);
    } // while

    return EXIT_SUCCESS;
} // overlap_main


typedef struct
{
    gva_uint start;
    gva_uint end;
    gva_uint inserted;
    gva_uint distance;
    gva_uint label;
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

    Trie labels = trie_init(gva_std_allocator);
    Trie sequences = trie_init(gva_std_allocator);
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
        } // if

        ARRAY_APPEND(gva_std_allocator, entries,
            ((Entry)
            {
                .start = variant.start,
                .end = variant.end,
                .inserted = trie_insert(&sequences, variant.sequence.len, variant.sequence.str),
                .distance = distance,
                .label = trie_insert(&labels, id.len, id.str),
            }));
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

            GVA_LCS_Graph lhs_graph = {NULL};
            size_t* lhs_dels = NULL;
            size_t* lhs_as = NULL;
            size_t* lhs_cs = NULL;
            size_t* lhs_gs = NULL;
            size_t* lhs_ts = NULL;

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

                size_t const start = MIN(lhs.start, rhs.start);
                size_t const end = MAX(lhs.end, rhs.end);

                size_t const len_lhs = (lhs.start - start) + lhs.sequence.len + (end - lhs.end);
                size_t const len_rhs = (rhs.start - start) + rhs.sequence.len + (end - rhs.end);

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
                    GVA_String observed_lhs = gva_string_init(gva_std_allocator, len_lhs);
                    GVA_String observed_rhs = gva_string_init(gva_std_allocator, len_rhs);
                    if (observed_lhs.str == NULL || observed_rhs.str == NULL)
                    {
                        return EXIT_FAILURE;  // FIXME: OOM
                    } // if

                    memcpy((char*) observed_lhs.str, reference.str + start, lhs.start - start);
                    memcpy((char*) observed_lhs.str + lhs.start - start, lhs.sequence.str, lhs.sequence.len);
                    memcpy((char*) observed_lhs.str + lhs.start - start + lhs.sequence.len, reference.str + lhs.end, end - lhs.end);

                    memcpy((char*) observed_rhs.str, reference.str + start, rhs.start - start);
                    memcpy((char*) observed_rhs.str + rhs.start - start, rhs.sequence.str, rhs.sequence.len);
                    memcpy((char*) observed_rhs.str + rhs.start - start + rhs.sequence.len, reference.str + rhs.end, end - rhs.end);

                    distance = gva_edit_distance(gva_std_allocator, observed_lhs.len, observed_lhs.str, observed_rhs.len, observed_rhs.str);
                    gva_string_destroy(gva_std_allocator, observed_rhs);
                    gva_string_destroy(gva_std_allocator, observed_lhs);
                } // else

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

                if (lhs_graph.nodes == NULL)
                {
                    lhs_graph = gva_lcs_graph_init(gva_std_allocator, lhs.end - lhs.start, reference.str + lhs.start, lhs.sequence.len, lhs.sequence.str, lhs.start);
                    size_t const len = end - start + 1;
                    lhs_dels = bitset_init(gva_std_allocator, len);
                    lhs_as = bitset_init(gva_std_allocator, len);
                    lhs_cs = bitset_init(gva_std_allocator, len);
                    lhs_gs = bitset_init(gva_std_allocator, len);
                    lhs_ts = bitset_init(gva_std_allocator, len);
                    gva_lcs_graph_uniq_atomics(lhs_graph, lhs.start, lhs.start, lhs.end, lhs_dels, lhs_as, lhs_cs, lhs_gs, lhs_ts);
                } // if

                GVA_LCS_Graph rhs_graph = gva_lcs_graph_init(gva_std_allocator, rhs.end - rhs.start, reference.str + rhs.start, rhs.sequence.len, rhs.sequence.str, rhs.start);

                size_t const len = end - start + 1;
                size_t* rhs_dels = bitset_init(gva_std_allocator, len);
                size_t* rhs_as = bitset_init(gva_std_allocator, len);
                size_t* rhs_cs = bitset_init(gva_std_allocator, len);
                size_t* rhs_gs = bitset_init(gva_std_allocator, len);
                size_t* rhs_ts = bitset_init(gva_std_allocator, len);

                gva_lcs_graph_uniq_atomics(rhs_graph, lhs.start, rhs.start, rhs.end, rhs_dels, rhs_as, rhs_cs, rhs_gs, rhs_ts);

                bool const overlap = bitset_intersection_cnt(lhs_dels, rhs_dels) > 0 ||
                    bitset_intersection_cnt(lhs_as, rhs_as) > 0 ||
                    bitset_intersection_cnt(lhs_cs, rhs_cs) > 0 ||
                    bitset_intersection_cnt(lhs_gs, rhs_gs) > 0 ||
                    bitset_intersection_cnt(lhs_ts, rhs_ts) > 0;

                rhs_ts = bitset_destroy(gva_std_allocator, rhs_ts);
                rhs_gs = bitset_destroy(gva_std_allocator, rhs_gs);
                rhs_cs = bitset_destroy(gva_std_allocator, rhs_cs);
                rhs_as = bitset_destroy(gva_std_allocator, rhs_as);
                rhs_dels = bitset_destroy(gva_std_allocator, rhs_dels);

                fprintf(stderr, GVA_STRING_FMT " vs " GVA_STRING_FMT "\n",
                    GVA_STRING_PRINT(trie_string(labels, entries[i].label)),
                    GVA_STRING_PRINT(trie_string(labels, entries[j].label)));
                bool const disjoint = is_disjoint(gva_std_allocator, lhs_graph, rhs_graph);

                if (disjoint == overlap)
                {
                    fprintf(stderr, "ERROR: %d\n", disjoint);
                    fprintf(stderr, GVA_STRING_FMT " overlap " GVA_STRING_FMT " %zu\n",
                        GVA_STRING_PRINT(trie_string(labels, entries[i].label)),
                        GVA_STRING_PRINT(trie_string(labels, entries[j].label)),
                        (entries[i].distance + entries[j].distance - distance) / 2);
                        return EXIT_FAILURE;
                } // if

                gva_lcs_graph_destroy(gva_std_allocator, rhs_graph, false);

                if (overlap)
                {
                    fprintf(stdout, GVA_STRING_FMT " overlap " GVA_STRING_FMT " %zu\n",
                        GVA_STRING_PRINT(trie_string(labels, entries[i].label)),
                        GVA_STRING_PRINT(trie_string(labels, entries[j].label)),
                        (entries[i].distance + entries[j].distance - distance) / 2);
                } // if
                else
                {
                    fprintf(stdout, GVA_STRING_FMT " disjoint " GVA_STRING_FMT "\n",
                        GVA_STRING_PRINT(trie_string(labels, entries[i].label)),
                        GVA_STRING_PRINT(trie_string(labels, entries[j].label)));
                } // else
            } // for
            lhs_ts = bitset_destroy(gva_std_allocator, lhs_ts);
            lhs_gs = bitset_destroy(gva_std_allocator, lhs_gs);
            lhs_cs = bitset_destroy(gva_std_allocator, lhs_cs);
            lhs_as = bitset_destroy(gva_std_allocator, lhs_as);
            lhs_dels = bitset_destroy(gva_std_allocator, lhs_dels);
            gva_lcs_graph_destroy(gva_std_allocator, lhs_graph, false);
        } // for

        fprintf(stderr, "#combos: %zu\n", count);
    } // if

    entries = ARRAY_DESTROY(gva_std_allocator, entries);
    trie_destroy(&labels);
    trie_destroy(&sequences);

    gva_string_destroy(gva_std_allocator, reference);

    return EXIT_SUCCESS;
} // all_main


int
main(int argc, char* argv[static argc])
{
    // return allele_main(argc, argv);
    // return index_main(argc, argv);
    return overlap_main(argc, argv);
    // return all_main(argc, argv);
} // main
