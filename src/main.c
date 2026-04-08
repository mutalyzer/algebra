#include <errno.h>      // errno
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // FILE, stderr, stdout, fclose, fopen, fprintf
#include <stdlib.h>     // EXIT_*, atoll, qsort
#include <string.h>     // strerror, strlen

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
#include "common.h"         // MAX, MIN
#include "priority_queue.h" // Priority_Queue, priority_queue_*


#include <assert.h>


#define LINE_SIZE 8194

// #define REFERENCE_ID "NC_000022.11"
#define REFERENCE_ID "NC_000006.12"
// #define REFERENCE_ID "NC_000001.11"


static void
pq_dot_traverse(Priority_Queue const self, size_t const i, size_t const rhs_length)
{
    size_t const len = array_length(self.heap);
    size_t const lhs_idx = self.heap[i] / rhs_length;
    size_t const rhs_idx = self.heap[i] % rhs_length;
    fprintf(stderr, "%zu[label=\"{%zu, %zu}\\n%u %u\\n%u\"]\n",
        i, lhs_idx, rhs_idx,
        self.states[self.heap[i]].included, self.states[self.heap[i]].excluded,
        self.heap[i]);
    if (2 * i + 1 < len)
    {
        fprintf(stderr, "%zu->%zu\n", i, 2 * i + 1);
        pq_dot_traverse(self, 2 * i + 1, rhs_length);
        if (2 * i + 2 < len)
        {
            fprintf(stderr, "%zu->%zu\n", i, 2 * i + 2);
            pq_dot_traverse(self, 2 * i + 2, rhs_length);
        } // if
    } // if
} // pq_dot_traverse


static inline void
pq_dot(Priority_Queue const self, size_t const rhs_length)
{
    fprintf(stderr, "strict digraph{\nnode[fixedsize=true,shape=circle,width=1]\n");
    if (array_length(self.heap) > 0)
    {
        pq_dot_traverse(self, 0, rhs_length);
    } // if
    fprintf(stderr, "}\n");
} // pq_dot


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


// DEBUG / check
int
extract_main(int argc, char* argv[static argc])
{
    (void) argv;

    size_t line_count = 0;
    static char line[LINE_SIZE] = {0};
    while (fgets(line, sizeof(line), stdin) != NULL)
    {
        line_count += 1;

        size_t const len_ref = strcspn(line, "\t ");
        GVA_String const reference = {len_ref, line};
        size_t const len_obs = strcspn(line + len_ref + 1, "\n\t ");
        GVA_String const observed = {len_obs, line + len_ref + 1};

        fprintf(stdout, "%zu: " GVA_STRING_FMT " (%zu) " GVA_STRING_FMT " (%zu)\n", line_count,
            GVA_STRING_PRINT(reference), reference.len, GVA_STRING_PRINT(observed), observed.len);

        GVA_LCS_Graph graph = gva_lcs_graph_init(gva_std_allocator,
            reference.len, reference.str, observed.len, observed.str, 0);

        // fprintf(stdout, "    distance:   %zu\n", gva_lcs_graph_distance(graph));
        // fprintf(stdout, "    #nodes:     %zu\n", array_length(graph.nodes));
        // for (size_t i = 0; i < array_length(graph.nodes); ++i)
        // {
        //     fprintf(stdout, "        %zu: (%u, %u, %u)\n", i,
        //         graph.nodes[i].match.row, graph.nodes[i].match.col, graph.nodes[i].match.length);
        // } // for
        // fprintf(stdout, "    #edges:     %zu\n", array_length(graph.edges));
        // fprintf(stdout, "    #dom_nodes: %zu\n", array_length(graph.dom_nodes));
        // if (array_length(graph.dom_nodes) < 3)
        // {
        //     continue;
        // }

        // for (size_t i = 0; i < array_length(graph.dom_nodes); ++i)
        // {
        //     if (i > 0)
        //     {
        //         fprintf(stdout, "           " GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(gva_lcs_graph_local_supremal(graph, i - 1, i)));
        //     } // if
        //     fprintf(stdout, "        %zu: (%u, %u, %u) %u\n", i,
        //         graph.dom_nodes[i].match.row, graph.dom_nodes[i].match.col, graph.dom_nodes[i].match.length,
        //         graph.dom_nodes[i].distance);
        // } // for

        // fprintf(stdout, "    supremal:   " GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(gva_lcs_graph_supremal(graph)));

        gva_lcs_graph_dot(stdout, graph);
        GVA_Variant sup = gva_variant_dup(gva_std_allocator, gva_lcs_graph_supremal(graph));

        GVA_LCS_Graph graph2 = gva_lcs_graph_from_allele(gva_std_allocator, reference.len, reference.str, 1, &sup);
        gva_lcs_graph_dot(stdout, graph2);

        if (array_length(graph.dom_nodes) != array_length(graph2.dom_nodes))
        {
            fprintf(stderr, "length dom nodes mismatch\n");
            assert(0);
        } // if

        for (size_t i = 0; i < array_length(graph.dom_nodes); ++i)
        {
            if (graph.dom_nodes[i].match.row != graph2.dom_nodes[i].match.row ||
                graph.dom_nodes[i].match.col != graph2.dom_nodes[i].match.col ||
                graph.dom_nodes[i].match.length != graph2.dom_nodes[i].match.length)
            {
                fprintf(stderr, "dom nodes mismatch\n");
                assert(0);
            } // if
        } // for

        gva_string_destroy(gva_std_allocator, sup.sequence);
        gva_lcs_graph_destroy(gva_std_allocator, graph2, true);
        gva_lcs_graph_destroy(gva_std_allocator, graph, false);
    } // while

    return EXIT_SUCCESS;
} // extract_main


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


static inline size_t
variant_included(size_t const lhs_count, GVA_Variant const lhs,
    size_t const rhs_count, GVA_Variant const rhs)
{
    if (rhs.start > lhs.end + lhs_count - 1 ||
        lhs.start > rhs.end + rhs_count - 1)
    {
        return 0;
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

    if (lhs.sequence.len == 0 || rhs.sequence.len == 0)
    {
        return end - start;
    } // if
    if (lhs.sequence.len == 1 && rhs.sequence.len == 1)
    {
        return (lhs.sequence.str[0] == rhs.sequence.str[0]) + end - start;
    } // if

    size_t const lcs = (lhs.sequence.len + rhs.sequence.len -
        gva_edit_distance(gva_std_allocator, lhs.sequence.len, lhs.sequence.str + lhs_offset,
                                             rhs.sequence.len, rhs.sequence.str + rhs_offset)) / 2;
    return lcs + end - start;
} // variant_included


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

        static char filename[128] = {0};
        snprintf(filename, 128, "overlap/graph_%zu.dot", line_count);
        FILE* stream = fopen(filename, "w");
        gva_lcs_graph_dot(stream, lhs);
        gva_lcs_graph_dot(stream, rhs);
        fclose(stream);

        fprintf(stderr, "%sdistances lhs: %zu rhs: %zu :: ", line, gva_lcs_graph_distance(lhs), gva_lcs_graph_distance(rhs));
        {
            size_t const distance = gva_edit_distance(gva_std_allocator,
                gva_lcs_graph_supremal(lhs).sequence.len, gva_lcs_graph_supremal(lhs).sequence.str,
                gva_lcs_graph_supremal(rhs).sequence.len, gva_lcs_graph_supremal(rhs).sequence.str);
            fprintf(stderr, "%zu\n", distance);
        } // empty
        fprintf(stderr, "nodes lhs: %zu rhs %zu :: %zu\n", array_length(lhs.nodes), array_length(rhs.nodes), array_length(lhs.nodes) * array_length(rhs.nodes));

        size_t const distance = gva_lcs_graph_distance(lhs) + gva_lcs_graph_distance(rhs);
        Priority_Queue fringe = priority_queue_init(gva_std_allocator, array_length(lhs.nodes) * array_length(rhs.nodes));
        priority_queue_push(&fringe, lhs.source * array_length(rhs.nodes) + rhs.source, 0, 0);

        size_t steps = 0;
        while (!priority_queue_empty(fringe))
        {
            steps += 1;
            size_t const head = priority_queue_pop(&fringe);

            size_t const lhs_idx = head / array_length(rhs.nodes);
            size_t const rhs_idx = head % array_length(rhs.nodes);

            fprintf(stderr, "{%zu, %zu} :: %u %u  %zu\n", lhs_idx, rhs_idx, fringe.states[head].included, fringe.states[head].excluded, distance - 2 * fringe.states[head].included - fringe.states[head].excluded);
            if (distance - 2 * fringe.states[head].included - fringe.states[head].excluded == 0)
            {
                break;
            } // if

            if (lhs.nodes[lhs_idx].lambda != GVA_NULL)
            {
                priority_queue_push(&fringe, lhs.nodes[lhs_idx].lambda * array_length(rhs.nodes) + rhs_idx,
                    fringe.states[head].included, fringe.states[head].excluded);
            } // if

            if (rhs.nodes[rhs_idx].lambda != GVA_NULL)
            {
                priority_queue_push(&fringe, lhs_idx * array_length(rhs.nodes) + rhs.nodes[rhs_idx].lambda,
                    fringe.states[head].included, fringe.states[head].excluded);
            } // if

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

                    size_t const included = variant_included(lhs_count, lhs_variant, rhs_count, rhs_variant);
                    size_t const excluded = gva_variant_length(lhs_variant) + gva_variant_length(rhs_variant) - 2 * included;

/*
                    fprintf(stderr, "  {%u, %u} ", lhs.edges[i].tail, rhs.edges[j].tail);
                    fprintf(stderr, GVA_VARIANT_FMT " x %zu (%zu) vs " GVA_VARIANT_FMT " x %zu (%zu) ", GVA_VARIANT_PRINT(lhs_variant), lhs_count, gva_variant_length(lhs_variant), GVA_VARIANT_PRINT(rhs_variant), rhs_count, gva_variant_length(rhs_variant));
                    fprintf(stderr, "+ %zu %zu\n", included, excluded);
*/

                    // follow two edges
                    priority_queue_push(&fringe, lhs.edges[i].tail * array_length(rhs.nodes) + rhs.edges[j].tail,
                        fringe.states[head].included + included, fringe.states[head].excluded + excluded);

                    // follow only one edge
                    if (i == lhs.nodes[lhs_idx].edges)
                    {
                        priority_queue_push(&fringe, lhs_idx * array_length(rhs.nodes) + rhs.edges[j].tail,
                            fringe.states[head].included, fringe.states[head].excluded + gva_variant_length(rhs_variant));
                    } // if
                } // for

                // follow only one edge
                priority_queue_push(&fringe, lhs.edges[i].tail * array_length(rhs.nodes) + rhs_idx,
                    fringe.states[head].included, fringe.states[head].excluded + gva_variant_length(lhs_variant));
            } // for
/*
            pq_dot(fringe, array_length(rhs.nodes));
            for (size_t i = 0; i < fringe.capacity; ++i)
            {
                if (fringe.states[i].idx != GVA_NULL)
                {
                    size_t const lhs_idx = i / array_length(rhs.nodes);
                    size_t const rhs_idx = i % array_length(rhs.nodes);
                    fprintf(stderr, "%4zu: {%zu, %zu} :: %2u %2u %2d\n", i, lhs_idx, rhs_idx, fringe.states[i].included, fringe.states[i].excluded, fringe.states[i].idx);
                } // if
            } // for
*/
        } // while
        fprintf(stderr, "steps: %zu\nin queue: %zu\n", steps, array_length(fringe.heap));

        priority_queue_destroy(&fringe);

        gva_lcs_graph_destroy(gva_std_allocator, rhs, true);
        gva_lcs_graph_destroy(gva_std_allocator, lhs, true);
    } // while

    return EXIT_SUCCESS;
} // overlap_main


int
main(int argc, char* argv[static argc])
{
    // return allele_main(argc, argv);
    // return extract_main(argc, argv);
    // return index_main(argc, argv);
    return overlap_main(argc, argv);
    // return dot_main(argc, argv);
} // main
