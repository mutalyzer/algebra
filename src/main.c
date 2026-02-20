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
#include "../include/index.h"       // GVA_Index, gva_index_*, GVA_Result
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

// #define REFERENCE_ID "NC_000006.12"
#define REFERENCE_ID "NC_000001.11"


inline void
serialize_array(FILE* restrict const stream,
    void* restrict const self, size_t const item_size)
{
    size_t const length = array_length(self);
    if (fwrite(&length, sizeof(length), 1, stream) != 1 ||
        fwrite(self, item_size, length, stream) != length)
    {
        fprintf(stderr, "error: fwrite()\n");
        return;
    } // if
} // serialize_array


inline void*
deserialize_array(GVA_Allocator const allocator, FILE* const stream,
    size_t const item_size)
{
    size_t length = 0;
    if (fread(&length, sizeof(length), 1, stream) != 1)
    {
        fprintf(stderr, "error: fread()\n");
        return NULL;
    } // if

    if (length == 0)
    {
        return NULL;
    } // if

    void* const restrict array = array_init(allocator, length, item_size);
    if (array == NULL)
    {
        fprintf(stderr, "error: OOM\n");
        return NULL;
    } // if

    if (fread(array, item_size, length, stream) != length)
    {
        fprintf(stderr, "error: fread()\n");
        return allocator.allocate(allocator.context, array_header(array), length * item_size, 0);
    } // if

    array_header(array)->length = length;
    return array;
} // deserialize_array


void
lcs_graph_dot(FILE* const stream, GVA_LCS_Graph const graph)
{
    fprintf(stream, "strict digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\ni[label=\"\",shape=none,width=0]\ni->%u\n", graph.source);
    for (size_t i = 0; i < array_length(graph.nodes); ++i)
    {
        fprintf(stream, "%zu[label=\"(%u, %u, %u)\"%s]\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, graph.nodes[i].edges == GVA_NULL ? ",peripheries=2" : "");
        if (graph.nodes[i].lambda != GVA_NULL)
        {
            fprintf(stream, "%zu->%u[label=\"&lambda;\",style=dashed]\n", i, graph.nodes[i].lambda);
        } // if
        for (gva_uint j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
        {
            GVA_Variant variant;
            gva_uint const count = gva_edges(graph.observed.str,
                graph.nodes[i], graph.nodes[graph.edges[j].tail],
                i == graph.source, graph.nodes[graph.edges[j].tail].edges == GVA_NULL,
                &variant);
            if (count > 1)
            {
                fprintf(stream, "%zu->%u[label=\"" GVA_VARIANT_FMT " x %u\",penwidth=2]\n", i, graph.edges[j].tail, GVA_VARIANT_PRINT(variant), count);
            } // if
            else
            {
                fprintf(stream, "%zu->%u[label=\"" GVA_VARIANT_FMT "\"]\n", i, graph.edges[j].tail, GVA_VARIANT_PRINT(variant));
            } // else
        } // for
    } // for
    for (size_t i = 0; i < array_length(graph.dom_nodes); ++i)
    {
        fprintf(stream, "%u[penwidth=2]\n", graph.dom_nodes[i].link);
    } // for
    fprintf(stream, "}\n");
} // lcs_graph_dot


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

    errno = 0;
    if (fread(&reference.len, sizeof(reference.len), 1, stream) != 1)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return (GVA_String) {0, NULL};
    }  // if

    reference.str = allocator.allocate(allocator.context, NULL, 0, reference.len);
    if (reference.str == NULL)
    {
        fprintf(stderr, "OOM\n");
        return (GVA_String) {0, NULL};
    } // if

    errno = 0;
    if (fread((char*) reference.str, 1, reference.len, stream) != reference.len)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        gva_string_destroy(allocator, reference);
        return (GVA_String) {0, NULL};
    } // if

    return reference;
} // gva_fasta_sequence_blob


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
        rhs_obs = allocator.allocate(allocator.context, rhs_obs, rhs_len, 0);
        lhs_obs = allocator.allocate(allocator.context, lhs_obs, lhs_len, 0);
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

    size_t const rhs_distance = graph.dom_nodes[part_idx + 1].distance;
    if (lhs_distance >= rhs_distance)
    {
        return 1;  // overlap
    } // if

    GVA_Variant const lhs = construct_variant(allocator, len_ref, reference, tree, trie, nodes);

    GVA_Variant rhs;
    gva_edges(graph.observed.str,
              graph.dom_nodes[part_idx], graph.dom_nodes[part_idx + 1],
              part_idx == 0, part_idx == array_length(graph.dom_nodes) - 2,
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


int vcf_main(int argc, char* argv[static argc])
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
            GVA_LCS_Graph graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, array_length(variants) - last, variants + last);

            for (size_t i = 0; i < array_length(graph.dom_nodes) - 1; ++i)
            {
                GVA_Variant variant;
                gva_edges(graph.observed.str,
                      graph.dom_nodes[i], graph.dom_nodes[i + 1],
                      i == 0, i == array_length(graph.dom_nodes) - 2, &variant);
                //fprintf(stdout, "%u " GVA_VARIANT_FMT_SPDI "\n", graph.dom_nodes[i + 1].distance, GVA_VARIANT_PRINT_SPDI("NC_000001.11", variant));
                distance += graph.dom_nodes[i + 1].distance;
                ARRAY_APPEND(gva_std_allocator, lss, gva_variant_dup(gva_std_allocator, variant));
            } // for

            last = array_length(variants);

            gva_string_destroy(gva_std_allocator, graph.observed);
            gva_lcs_graph_destroy(gva_std_allocator, graph);
        } // if

        ARRAY_APPEND(gva_std_allocator, variants, gva_variant_dup(gva_std_allocator, trimmed));
        line_count += 1;
    } // while
    if (last < array_length(variants))
    {
        GVA_LCS_Graph graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, array_length(variants) - last, variants + last);

        for (size_t i = 0; i < array_length(graph.dom_nodes) - 1; ++i)
        {
            GVA_Variant variant;
            gva_edges(graph.observed.str,
                  graph.dom_nodes[i], graph.dom_nodes[i + 1],
                  i == 0, i == array_length(graph.dom_nodes) - 2, &variant);
            //fprintf(stdout, "%u " GVA_VARIANT_FMT_SPDI "\n", graph.dom_nodes[i + 1].distance, GVA_VARIANT_PRINT_SPDI("NC_000001.11", variant));
            distance += graph.dom_nodes[i + 1].distance;
            ARRAY_APPEND(gva_std_allocator, lss, gva_variant_dup(gva_std_allocator, variant));
        } // for

        last = array_length(variants);

        gva_string_destroy(gva_std_allocator, graph.observed);
        gva_lcs_graph_destroy(gva_std_allocator, graph);
    } // if

    fprintf(stderr, "===INPUT===\n");
    fprintf(stderr, "#variants: %zu\n", array_length(variants));
    fprintf(stderr, "#dropped:  %zu\n", dropped);

    GVA_String observed = gva_patch(gva_std_allocator, reference.len, reference.str, array_length(variants), variants);

    fasta_blob_write(stdout, observed);

    return EXIT_SUCCESS;
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
    *nodes = ARRAY_DESTROY(gva_std_allocator, *nodes);  //  FIXME
} // repair_is_contained


static bool
parse_line(char const line[static LINE_SIZE],
    size_t* id_len, GVA_Variant* variant, size_t* distance)
{
    *id_len = strcspn(line, "\t ");
    if (*id_len == 0)
    {
        return false;
    } // if

    size_t const spdi_len = strcspn(line + *id_len + 1, "\t ");
    if (gva_parse_spdi(spdi_len, line + *id_len + 1, variant) == 0)
    {
        return false;
    } // if

    *distance = parse_number(line, &(size_t) {*id_len + spdi_len + 2});
    if (*distance == 0)
    {
        return false;
    } // if

    return true;
} // parse_line


static bool
parse_line2(char const line[static LINE_SIZE],
    size_t* id_len, GVA_Variant* variant)
{
    *id_len = strcspn(line, "\t ");
    if (*id_len == 0)
    {
        fprintf(stderr, "here?!\n");
        return false;
    } // if

    size_t const spdi_len = strcspn(line + *id_len + 1, "\n");
    if (gva_parse_spdi(spdi_len, line + *id_len + 1, variant) == 0)
    {
        fprintf(stderr, "or here?!\n");
        return false;
    } // if

    return true;
} // parse_line2
//

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

    GVA_String reference = {0, NULL};
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
        if (!parse_line2(line, &len_id, &variant))
        {
            fprintf(stderr, "parsing failed at line %zu: %s\n", line_count, line);
            continue;
        } // if

        GVA_LCS_Graph graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, 1, &variant);

        // fprintf(stderr, "\nQuery (" GVA_STRING_FMT "): " GVA_VARIANT_FMT_SPDI " (%u)\n", GVA_STRING_PRINT(((GVA_String) {id_len, line})), GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, graph.supremal), graph.distance);
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
        gva_string_destroy(gva_std_allocator, graph.observed);
        gva_lcs_graph_destroy(gva_std_allocator, graph);
    } // while

    index = gva_index_destroy(index);
    gva_string_destroy(gva_std_allocator, reference);

    return EXIT_SUCCESS;
} // index_main


// FIXME: to string module?
static char*
strrev(size_t const n, char str[static n])
{
    size_t i = n - 1;
    size_t j = 0;
    while (i > j)
    {
        char const ch = str[i];
        str[i] = str[j];
        str[j] = ch;
        i -= 1;
        j += 1;
    } // while
    return str;
} // strrev


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

        GVA_Variant trimmed = prefix_trimmed(reference.len, reference.str, variant);

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

    reference = strrev(len_ref, (char*) reference);
    observed = strrev(len_obs, (char*) observed);

    GVA_Matches backward = gva_edit_distance_matches(gva_std_allocator, len_ref, reference, len_obs, observed);

    if (forward.distance != backward.distance)
    {
         fprintf(stderr, "distance mismatch between f and b\n");
         exit(EXIT_FAILURE);  // FIXME
    } // if

    reference = strrev(len_ref, (char*) reference);
    observed = strrev(len_obs, (char*) observed);

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
                graph.dom_nodes[i], graph.dom_nodes[i + 1],
                i == 0, i == array_length(graph.dom_nodes) - 2,
                &variant);
            printf(GVA_VARIANT_FMT_SPDI " %u\n", GVA_VARIANT_PRINT_SPDI(REFERENCE_ID, variant), graph.dom_nodes[i + 1].distance);
        } // for
        gva_lcs_graph_destroy(gva_std_allocator, graph);
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
    // lcs_graph_dot(stderr, graph);
    for (size_t i = 0; i < array_length(graph.dom_nodes) - 1; ++i)
    {
        GVA_Variant part;
        gva_edges(graph.observed.str,
            graph.dom_nodes[i], graph.dom_nodes[i + 1],
            i == 0, i == array_length(graph.dom_nodes) - 2,
            &part);
        fprintf(stderr, GVA_VARIANT_FMT " %u\n", GVA_VARIANT_PRINT(part), graph.dom_nodes[i + 1].distance);
    } // for
    fprintf(stderr, GVA_VARIANT_FMT " %u\n", GVA_VARIANT_PRINT(graph.supremal), graph.distance);

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
    errno = 0;
    FILE *stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    GVA_String reference = {0, NULL};
    reference = gva_fasta_sequence(gva_std_allocator, stream);
    fclose(stream);

    fprintf(stderr, "reference length: %zu\n", reference.len);

    GVA_LCS_Graph graph = gva_lcs_graph_init(gva_std_allocator, 31, "CAGATGACAGGGTTGGGCTCAGAGTCAGAGT", 25, "CTCTCATACACACTCTCATCAGATT", 0);
    fprintf(stderr, "graph distance: %u\n", graph.distance);

    lcs_graph_dot(stderr, graph);

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
    FILE *stream = fopen(argv[1], "r");
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
    FILE *stream = fopen(argv[1], "r");
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
    // return slice_blob_main(argc, argv);
    // return fasta_blob_write(argc, argv);
    // return vcf_main(argc, argv);
    // return dbsnp_main(argc, argv);
    // return locals_main(argc, argv);
    // return extract_main(argc, argv);
    // return make_ref_blob_main(argc, argv);
    // return make_obs_blob_main(argc, argv);
    return index_main(argc, argv);
} // main
