#include <errno.h>      // errno
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // FILE, stderr, stdout, fclose, fopen, fprintf
#include <stdlib.h>     // EXIT_*, atoll, qsort
#include <string.h>     // strerror, strlen

#include "../include/index.h"       // GVA_Index, gva_index_*, GVA_Query_Result
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_lcs_graph_*
#include "../include/relations.h"   // GVA_RELATION_LABELS
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/string.h"      // GVA_String, gva_string_destroy
#include "../include/utils.h"       // gva_fasta_sequence*, gva_lcs_graph_dot
#include "../include/variant.h"     // GVA_VARIANT_*, GVA_Variant, gva_parse_spdi, gva_variant_*
#include "array.h"          // ARRAY_*, array_length
#include "common.h"         // MAX, MIN


#define LINE_SIZE 8194

// #define REFERENCE_ID "NC_000022.11"
// #define REFERENCE_ID "NC_000006.12"
#define REFERENCE_ID "NC_000001.11"


// line: alphanumeric_id SPDI [distance]
static bool
parse_line(char const line[static LINE_SIZE],
    GVA_String* const id, GVA_Variant* const variant, size_t* const distance)
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
    struct GVA_Query_Allele const lhs = *(struct GVA_Query_Allele*) a;
    struct GVA_Query_Allele const rhs = *(struct GVA_Query_Allele*) b;

    if (lhs.included > rhs.included ||
        (lhs.included == rhs.included && lhs.excluded < rhs.excluded))
    {
        return -1;
    } // if
    return 1;
} // compare_query_alleles


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
            fprintf(stderr, "parsing failed at line %zu: %s\n", line_count, line);
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
            fprintf(stderr, "parsing failed at line %zu: %s\n", line_count, line);
            continue;
        } // if
        GVA_LCS_Graph graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, 1, &variant);

        fprintf(stderr, "\nQuery (" GVA_STRING_FMT "): " GVA_VARIANT_FMT_SPDI " (%zu)\n", GVA_STRING_PRINT(id), GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, gva_lcs_graph_supremal(graph)), gva_lcs_graph_distance(graph));
        GVA_Query_Result result = gva_index_query(gva_std_allocator, index, graph);

        if (array_length(result.alleles) > 0)
        {
            qsort(result.alleles, array_length(result.alleles), sizeof(*result.alleles), compare_query_alleles);
        } // if

        for (size_t i = 0; i < array_length(result.alleles); ++i)
        {
            fprintf(stdout, GVA_STRING_FMT " %s " GVA_STRING_FMT " %u %u %zu\n",
                GVA_STRING_PRINT(gva_index_id(index, result.alleles[i].idx)),
                GVA_RELATION_LABELS[result.alleles[i].relation],
                GVA_STRING_PRINT(id),
                result.alleles[i].included, result.alleles[i].excluded, array_length(graph.dom_nodes) - 1);
            fprintf(stderr, "    " GVA_STRING_FMT ": %u %u %s\n",
                GVA_STRING_PRINT(gva_index_id(index, result.alleles[i].idx)),
                result.alleles[i].included, result.alleles[i].excluded,
                GVA_RELATION_LABELS[result.alleles[i].relation]);
            for (size_t j = result.alleles[i].hits.start; j < result.alleles[i].hits.end; ++j)
            {
                fprintf(stderr, "        [%u, %u) [%u, %u): %u %u %s\n",
                    result.hits[j].query.start, result.hits[j].query.end,
                    result.hits[j].index.start, result.hits[j].index.end,
                    result.hits[j].included, result.hits[j].excluded,
                    GVA_RELATION_LABELS[result.hits[j].relation]);
            } // for
        } // for
        result.alleles = ARRAY_DESTROY(gva_std_allocator, result.alleles);
        result.hits = ARRAY_DESTROY(gva_std_allocator, result.hits);
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

        fprintf(stdout, "    distance:   %zu\n", gva_lcs_graph_distance(graph));
        fprintf(stdout, "    #nodes:     %zu\n", array_length(graph.nodes));
        for (size_t i = 0; i < array_length(graph.nodes); ++i)
        {
            fprintf(stdout, "        %zu: (%u, %u, %u)\n", i,
                graph.nodes[i].match.row, graph.nodes[i].match.col, graph.nodes[i].match.length);
        } // for
        fprintf(stdout, "    #edges:     %zu\n", array_length(graph.edges));
        fprintf(stdout, "    #dom_nodes: %zu\n", array_length(graph.dom_nodes));
        for (size_t i = 0; i < array_length(graph.dom_nodes); ++i)
        {
            if (i > 0)
            {
                fprintf(stdout, "           " GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(gva_lcs_graph_local_supremal(graph, i - 1, i)));
            } // if
            fprintf(stdout, "        %zu: (%u, %u, %u) %u\n", i,
                graph.dom_nodes[i].match.row, graph.dom_nodes[i].match.col, graph.dom_nodes[i].match.length,
                graph.dom_nodes[i].distance);
        } // for

        fprintf(stdout, "    supremal:   " GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(gva_lcs_graph_supremal(graph)));

        gva_lcs_graph_dot(stdout, graph);

        gva_lcs_graph_destroy(gva_std_allocator, graph, false);
    } // while

    return EXIT_SUCCESS;
} // extract_main


int
main(int argc, char* argv[static argc])
{
    // return allele_main(argc, argv);
    // return extract_main(argc, argv);
    return index_main(argc, argv);
} // main
