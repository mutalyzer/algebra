#include <errno.h>      // errno
#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // FILE, stderr, stdout, fclose, fopen, fprintf
#include <stdlib.h>     // EXIT_*, atoll
#include <string.h>     // strerror, strlen

#include "../include/edit.h"        // gva_edit_distance
#include "../include/index.h"       // GVA_Index, gva_index_*, GVA_Result
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, GVA_Variant, gva_lcs_graph_*, gva_edges
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/string.h"      // GVA_String, gva_string_destroy
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "../include/utils.h"       // gva_fasta_sequence*, gva_lcs_graph_dot
#include "../include/variant.h"     // GVA_VARIANT_*, GVA_Variant, gva_parse_spdi
#include "array.h"          // ARRAY_*, array_length
#include "bitset.h"         // bitset_*
#include "common.h"         // MAX, MIN
#include "hash_table.h"     // GVA_NOT_FOUND, HASH_TABLE_KEY, hash_table_*


#define LINE_SIZE 8194

// #define REFERENCE_ID "NC_000006.12"
#define REFERENCE_ID "NC_000001.11"


static void
fasta_blob_write(FILE* stream, GVA_String const sequence)
{
    errno = 0;
    if (fwrite(&sequence.len, sizeof(sequence.len), 1, stream) != 1)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return;
    } // if

    errno = 0;
    if (fwrite(sequence.str, 1, sequence.len, stream) != sequence.len)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return;
    } // if
} // fasta_blob_write


// line: alphanumeric_id SPDI [distance]
static bool
parse_line(char const line[static LINE_SIZE],
    size_t* const len_id, GVA_Variant* const variant, size_t* const distance)
{
    *len_id = strcspn(line, "\t ");
    if (*len_id == 0)
    {
        return false;
    } // if

    size_t const len_spdi = strcspn(line + *len_id + 1, "\n\t ");
    if (gva_parse_spdi(len_spdi, line + *len_id + 1, variant) == 0)
    {
        return false;
    } // if

    *distance = atoll(line + *len_id + len_spdi + 2);
    return true;
} // parse_line


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
        size_t len_id = 0;
        GVA_Variant variant;
        size_t distance = 0;
        if (!parse_line(line, &len_id, &variant, &distance))
        {
            fprintf(stderr, "parsing failed at line %zu: %s\n", line_count, line);
            continue;
        } // if
        gva_index_insert(index, len_id, line, variant, distance);
    } // while

    fclose(stream);

    line_count = 0;
    while (fgets(line, sizeof(line), stdin) != NULL)
    {
        line_count += 1;
        size_t len_id = 0;
        GVA_Variant variant;
        if (!parse_line(line, &len_id, &variant, &(size_t) {0}))
        {
            fprintf(stderr, "parsing failed at line %zu: %s\n", line_count, line);
            continue;
        } // if

        GVA_LCS_Graph graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, 1, &variant);

        fprintf(stderr, "\nQuery (" GVA_STRING_FMT "): " GVA_VARIANT_FMT_SPDI " (%zu)\n", GVA_STRING_PRINT(((GVA_String) {len_id, line})), GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, gva_lcs_graph_supremal(graph)), gva_lcs_graph_distance(graph));
        GVA_Query_Result* results = gva_index_query(gva_std_allocator, index, graph);
        for (size_t i = 0; i < array_length(results); ++i)
        {
            fprintf(stdout, GVA_STRING_FMT " %s " GVA_STRING_FMT " %zu %zu\n",
                GVA_STRING_PRINT(results[i].allele),
                GVA_RELATION_LABELS[results[i].relation],
                GVA_STRING_PRINT(((GVA_String) {len_id, line})),
                results[i].included, results[i].excluded);
        } // for

        ARRAY_DESTROY(gva_std_allocator, results);
        gva_lcs_graph_destroy(gva_std_allocator, graph, true);
    } // while

    index = gva_index_destroy(index);
    gva_string_destroy(gva_std_allocator, reference);

    return EXIT_SUCCESS;
} // index_main


GVA_String
vcf2obs(GVA_Allocator const allocator, GVA_String const reference, FILE* stream)
{
    GVA_Variant* variants = NULL;

    size_t line_count = 0;
    size_t dropped = 0;
    static char line[LINE_SIZE] = {0};
    while (fgets(line, sizeof(line), stream) != NULL)
    {
        size_t const len = strlen(line) - 1;
        GVA_Variant variant;
        if (gva_parse_spdi(len, line, &variant) == 0)
        {
            fprintf(stderr, "error: SPDI parsing failed at line %zu: %s", line_count + 1, line);
            continue;
        } // if

        GVA_Variant const trimmed = gva_variant_prefix_trimmed(reference.len, reference.str, variant);

        if (array_length(variants) > 0 && trimmed.start < variants[array_length(variants) - 1].end)
        {
            //fprintf(stderr, "dropped: at line %zu: %s", line_count + 1, line);
            dropped += 1;
            continue;
        } // if

        ARRAY_APPEND(allocator, variants, gva_variant_dup(allocator, trimmed));
        line_count += 1;
    } // while

    fprintf(stderr, "===INPUT===\n");
    fprintf(stderr, "#variants: %zu\n", array_length(variants));
    fprintf(stderr, "#dropped:  %zu\n", dropped);

    return gva_patch(allocator, reference.len, reference.str, array_length(variants), variants);
} // vcf2obs


//#define SMALL


void
local_supremal(size_t const len_ref, char const reference[static len_ref],
    size_t const len_obs, char const observed[static len_obs],
    size_t const offset)
{
    //fprintf(stderr, "LOCAL %zu %zu :: %zu\n", len_ref, len_obs, offset);
    if (len_ref == 0 || len_obs == 0)
    {
        // fprintf(stderr, "triv distance: %zu\n", len_ref + len_obs);
        printf(GVA_VARIANT_FMT_SPDI " %zu\n", GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, ((GVA_Variant) {offset, offset + len_ref, {len_obs, observed}})), len_ref + len_obs);
        return;
    } // if
    
    GVA_Matches forward = gva_edit_distance_matches(gva_std_allocator, len_ref, reference, len_obs, observed);

    gva_string_reverse((GVA_String) {len_ref, reference});
    gva_string_reverse((GVA_String) {len_obs, observed});

    GVA_Matches backward = gva_edit_distance_matches(gva_std_allocator, len_ref, reference, len_obs, observed);

    if (forward.distance != backward.distance)
    {
         fprintf(stderr, "distance mismatch between f and b\n");
         exit(EXIT_FAILURE);  // FIXME
    } // if

    gva_string_reverse((GVA_String) {len_ref, reference});
    gva_string_reverse((GVA_String) {len_obs, observed});

    size_t sum = 0;
    size_t prev_lcs_pos = 0;
    size_t prev_row = -1;
    size_t prev_col = -1;
    for (size_t i = 0; i < forward.max_lcs_pos; ++i)
    {
        size_t const j = forward.max_lcs_pos - i - 1;
        if (forward.matches[i].row == len_ref - backward.matches[j].row - 1 &&
            forward.matches[i].col == len_obs - backward.matches[j].col - 1 &&
            (forward.uniq[i] == 1 || backward.uniq[j] == 1))
        {
            #ifdef SMALL
            fprintf(stderr, "common: %zu: (%u, %u)\n", i, forward.matches[i].row, forward.matches[i].col);
            #endif
            if (i > prev_lcs_pos + 1 ||
                forward.matches[i].row > prev_row + 1 ||
                forward.matches[i].col > prev_col + 1)
            {
                size_t const distance = forward.matches[i].row + forward.matches[i].col - 2 * i - sum;
                if (distance > 0)
                {
                    sum += distance;
                    local_supremal(forward.matches[i].row - prev_row - 1, reference + prev_row + 1, forward.matches[i].col - prev_col - 1, observed + prev_col + 1, offset + prev_row + 1);
                } // if
            } // if
            prev_lcs_pos = i;
            prev_row = forward.matches[i].row;
            prev_col = forward.matches[i].col;
        } // if
    } // for
    if (len_ref > prev_row + 1 || len_obs > prev_col + 1)
    {
        size_t const distance = len_ref + len_obs - 2 * forward.max_lcs_pos - sum;
        sum += distance;
        GVA_LCS_Graph graph = gva_lcs_graph_init(gva_std_allocator, len_ref - prev_row - 1, reference + prev_row + 1, len_obs - prev_col - 1, observed + prev_col + 1, offset + prev_row + 1);
        for (size_t i = 0; i < array_length(graph.dom_nodes) - 1; ++i)
        {
            GVA_Variant variant;
            gva_edges(graph.observed.str,
                graph.dom_nodes[i].match, graph.dom_nodes[i + 1].match,
                i == 0, i + 1 == array_length(graph.dom_nodes) - 1,
                &variant);
            printf(GVA_VARIANT_FMT_SPDI " %u\n", GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, variant), graph.dom_nodes[i + 1].distance - graph.dom_nodes[i].distance);
        } // for
        gva_lcs_graph_destroy(gva_std_allocator, graph, false);
    } // if

    if (forward.distance != sum)
    {
         fprintf(stderr, "distance mismatch between f and sum\n");
         exit(EXIT_FAILURE);  // FIXME
    } // if

    backward.uniq = gva_std_allocator.allocate(gva_std_allocator.context, backward.uniq, MIN(len_ref, len_obs), 0);
    backward.matches = gva_std_allocator.allocate(gva_std_allocator.context, backward.matches, MIN(len_ref, len_obs) * sizeof(*backward.matches), 0);
    forward.uniq = gva_std_allocator.allocate(gva_std_allocator.context, forward.uniq, MIN(len_ref, len_obs), 0);
    forward.matches = gva_std_allocator.allocate(gva_std_allocator.context, forward.matches, MIN(len_ref, len_obs) * sizeof(*forward.matches), 0);
} // local_supremal


int
wu_main(int argc, char* argv[static argc])
{
    if (argc < 3)
    {
        fprintf(stderr, "usage: %s reference observed\n", argv[0]);
        return EXIT_FAILURE;
    } // if

#ifdef SMALL
    GVA_String reference = {strlen(argv[1]), argv[1]};
    GVA_String observed = {strlen(argv[2]), argv[2]};

    GVA_LCS_Graph graph = gva_lcs_graph_init(gva_std_allocator, reference.len, reference.str, observed.len, observed.str, 0);
    // gva_lcs_graph_dot(stderr, graph);
    for (size_t i = 0; i < array_length(graph.dom_nodes) - 1; ++i)
    {
        GVA_Variant part;
        gva_edges(graph.observed.str,
            graph.dom_nodes[i], graph.dom_nodes[i + 1],
            i == 0, i + 1 == array_length(graph.dom_nodes) - 1,
            &part);
        fprintf(stderr, GVA_VARIANT_FMT " %u\n", GVA_VARIANT_PRINT(part), graph.dom_nodes[i + 1].distance);
    } // for
    fprintf(stderr, GVA_VARIANT_FMT " %zu\n", GVA_VARIANT_PRINT(gva_lcs_graph_supremal(graph)), gva_lcs_graph_distance(graph));

    gva_lcs_graph_destroy(gva_std_allocator, graph);

    return 0;

#else
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

    // GVA_String observed = gva_fasta_sequence_blob(gva_std_allocator, stream);  // blob
    GVA_String observed = vcf2obs(gva_std_allocator, reference, stream);  // vcf
    fclose(stream);
    fprintf(stderr, "observed length: %zu\n", observed.len);

    if (reference.str == NULL || observed.str == NULL)
    {
        gva_string_destroy(gva_std_allocator, reference);
        gva_string_destroy(gva_std_allocator, observed);
        return EXIT_FAILURE;
    } // if
#endif

    local_supremal(reference.len, reference.str, observed.len, observed.str, 0);

#ifndef SMALL
    gva_string_destroy(gva_std_allocator, observed);
    gva_string_destroy(gva_std_allocator, reference);
#endif

    return EXIT_SUCCESS;
} // wu_main


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
                GVA_Variant variant = {0};
                gva_edges(graph.observed.str,
                    graph.dom_nodes[i - 1].match, graph.dom_nodes[i].match,
                    i - 1 == 0, i == array_length(graph.dom_nodes) - 1,
                    &variant);
                fprintf(stdout, "           " GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(variant));
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
make_ref_blob_main(int argc, char* argv[static argc])
{
    if (argc < 2)
    {
        fprintf(stderr, "usage: %s reference\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    errno = 0;
    FILE* stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    GVA_String reference = gva_fasta_sequence(gva_std_allocator, stream);
    fclose(stream);
    fprintf(stderr, "reference length: %zu\n", reference.len);

    fasta_blob_write(stdout, reference);

    return EXIT_SUCCESS;
} // make_ref_blob_main


int
make_obs_blob_main(int argc, char* argv[static argc])
{
    if (argc < 2)
    {
        fprintf(stderr, "usage: %s reference.blob\n", argv[0]);
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

    GVA_String observed = vcf2obs(gva_std_allocator, reference, stdin);
    fasta_blob_write(stdout, observed);

    return EXIT_SUCCESS;
} // make_obs_blob_main


int
main(int argc, char* argv[static argc])
{
    // return wu_main(argc, argv);
    // return fasta_blob_write(argc, argv);
    // return extract_main(argc, argv);
    // return make_ref_blob_main(argc, argv);
    // return make_obs_blob_main(argc, argv);
    return index_main(argc, argv);
} // main
