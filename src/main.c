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
#include "../include/stabbing.h"    // GVA_Stabbing_*, gva_stabbing_*
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/string.h"      // GVA_String, gva_string_destroy
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "../include/utils.h"       // gva_fasta_sequence
#include "../include/variant.h"     // GVA_VARIANT_*, GVA_Variant, gva_parse_spdi
#include "array.h"      // ARRAY_DESTROY, array_length
#include "bitset.h"     // bitset_*
#include "common.h"     // MAX, MIN
#include "trie.h"       // Trie, trie_*


#define TIC(x) clock_t const x = clock()
#define TOC(x) fprintf(stderr, "Elapsed time: %s (seconds): %.2f\n", (#x), (double)(clock() - (x)) / CLOCKS_PER_SEC)

#define LINE_SIZE 8194


static void
lcs_graph_svg(FILE* const stream, GVA_LCS_Graph const graph)
{
    gva_uint const offset = graph.nodes[graph.source].row;
    gva_uint const scale = 30;
    fprintf(stream, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 %u %u\">\n", (graph.nodes[0].col + graph.nodes[0].length) * scale + 2 * scale, (graph.nodes[0].row + graph.nodes[0].length - offset) * scale + 2 * scale);

    fprintf(stream, "<g stroke=\"black\">\n");
    for (size_t i = 0; i < array_length(graph.nodes); ++i)
    {
        if (graph.nodes[i].lambda != GVA_NULL)
        {
            gva_uint const row = (graph.nodes[i].row - offset) * scale + scale;
            gva_uint const col = graph.nodes[i].col * scale + scale;
            gva_uint const to_row = (graph.nodes[graph.nodes[i].lambda].row - offset) * scale + scale;
            gva_uint const to_col = graph.nodes[graph.nodes[i].lambda].col * scale + scale;
            fprintf(stream, "<line x1=\"%u\" y1=\"%u\" x2=\"%u\" y2=\"%u\" stroke-dasharray=\"2\"/>\n", col, row, to_col, to_row);
        } // if

        for (gva_uint j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
        {
            GVA_Variant variant;
            gva_uint const count = gva_edges(graph.observed.str,
                                             graph.nodes[i], graph.nodes[graph.edges[j].tail],
                                             i == graph.source, graph.nodes[graph.edges[j].tail].edges == GVA_NULL,
                                             &variant);

            gva_uint const row = (variant.start - (variant.start > 0) - offset); // * scale + scale;
            gva_uint const col = (variant.sequence.str - graph.observed.str - (variant.start > 0)); //  * scale + scale;
            gva_uint const to_row = (variant.end + count - 1 - offset); //  * scale + scale;
            gva_uint const to_col = (col + count - (variant.start == 0) + variant.sequence.len); // * scale + scale;

            fprintf(stream, "<line x1=\"%u\" y1=\"%u\" x2=\"%u\" y2=\"%u\" stroke-width=\"%.2f\"/>\n", col * scale + scale, row * scale + scale, to_col * scale + scale, to_row * scale + scale, (double) count / 2.0);
        } // for
    } // for
    fprintf(stream, "</g>\n");

    for (size_t i = 0; i < array_length(graph.nodes); ++i)
    {
        gva_uint const row = (graph.nodes[i].row - offset) * scale + scale;
        gva_uint const col = graph.nodes[i].col * scale + scale;
        if (graph.nodes[i].length < 2)
        {
            fprintf(stream, "<circle cx=\"%u\" cy=\"%u\" r=\"10\" fill=\"white\" stroke=\"black\" stroke-width=\".5\"/>\n", col, row);
        } // if
        else
        {
            gva_uint const to_row = (graph.nodes[i].row - offset + graph.nodes[i].length) * scale + scale;
            gva_uint const to_col = (graph.nodes[i].col + graph.nodes[i].length) * scale + scale;
            fprintf(stream, "<line x1=\"%u\" y1=\"%u\" x2=\"%u\" y2=\"%u\" stroke=\"black\" stroke-width=\"20.5\" stroke-linecap=\"round\"/>\n", col, row, to_col, to_row);
            fprintf(stream, "<line x1=\"%u\" y1=\"%u\" x2=\"%u\" y2=\"%u\" stroke=\"white\" stroke-width=\"19.5\" stroke-linecap=\"round\"/>\n", col, row, to_col, to_row);
        } // else
        fprintf(stream, "<text x=\"%u\" y=\"%u\" font-size=\"4px\" text-anchor=\"middle\">(%u, %u, %u)</text>\n", col, row + 1, graph.nodes[i].row - offset, graph.nodes[i].col, graph.nodes[i].length);
        if (graph.nodes[i].lambda != GVA_NULL)
        {
            gva_uint const row = (graph.nodes[graph.nodes[i].lambda].row - offset) * scale + scale;
            gva_uint const col = graph.nodes[graph.nodes[i].lambda].col * scale + scale;
            fprintf(stream, "<text x=\"%u\" y=\"%u\" font-size=\"4px\" text-anchor=\"middle\">(%u, %u, %u)</text>\n", col, row + 1, graph.nodes[graph.nodes[i].lambda].row - offset, graph.nodes[graph.nodes[i].lambda].col, graph.nodes[graph.nodes[i].lambda].length);
        } // ifr
    } // for
    fprintf(stream, "</svg>\n");
} // lcs_graph_svg


static void
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
    for (size_t i = 0; i < array_length(graph.local_supremal); ++i)
    {
        fprintf(stream, "%u[penwidth=2]\n", graph.local_supremal[i].lambda);
    } // for
    fprintf(stream, "}\n");
} // lcs_graph_dot


static void
bfs_traversal(FILE* const stream, GVA_LCS_Graph const graph)
{
    struct
    {
        gva_uint depth;
        gva_uint next;
    }* table = gva_std_allocator.allocate(gva_std_allocator.context, NULL, 0, sizeof(*table) * array_length(graph.nodes));
    if (table == NULL)
    {
        return;
    } // if

    for (size_t i = 0; i < array_length(graph.nodes); ++i)
    {
        table[i].depth = 0;
        table[i].next = GVA_NULL;
    } // for
    gva_uint head = graph.source;
    gva_uint tail = graph.source;

    bool first = true;
    fprintf(stream, "    \"edges\": [\n");
    while (head != GVA_NULL)
    {
        gva_uint len = graph.nodes[head].length;
        for (gva_uint i = head; i != GVA_NULL; i = graph.nodes[i].lambda)
        {
            if (i != head)
            {
                len += graph.nodes[i].length;
            } // if
            for (gva_uint j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
            {
                GVA_Variant variant;
                gva_uint const count = gva_edges(graph.observed.str,
                    graph.nodes[i], graph.nodes[graph.edges[j].tail],
                    i == graph.source, graph.nodes[graph.edges[j].tail].edges == GVA_NULL,
                    &variant);

                for (gva_uint k = 0; k < count; ++k)
                {
                    if (!first)
                    {
                        printf(",\n");
                    } // if
                    first = false;
                    fprintf(stream, "         {\"head\": \"s%u\", \"tail\": \"s%u\", \"variant\": \"" GVA_VARIANT_FMT "\"}", head, graph.edges[j].tail, GVA_VARIANT_PRINT(variant));
                    variant.start += 1;
                    variant.end += 1;
                    variant.sequence.str += 1;
                } // for

                if (table[graph.edges[j].tail].depth > 0)
                {
                    continue;
                } // if

                table[graph.edges[j].tail].depth = table[i].depth + 1;
                table[tail].next = graph.edges[j].tail;
                tail = graph.edges[j].tail;
            } // for
        } // for
        head = table[head].next;
    } // while

    fprintf(stream, "\n    ],\n");

    table = gva_std_allocator.allocate(gva_std_allocator.context, table, sizeof(*table) * array_length(graph.nodes), 0);
} // bfs_traversal


static void
lcs_graph_json(FILE* const stream, GVA_LCS_Graph const graph)
{
    fprintf(stream, "{\n    \"source\": \"s%u\",\n    \"nodes\": {\n", graph.source);
    for (gva_uint i = 0; i < array_length(graph.nodes); ++i)
    {
        fprintf(stream, "        \"s%u\": {\"row\": %u, \"col\": %u, \"length\": %u}%s\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, i < array_length(graph.nodes) - 1 ? "," : "");
    } // for
    fprintf(stream, "    },\n");
    bfs_traversal(stream, graph);
    fprintf(stream, "    \"local_supremal\": [\n");
    if (array_length(graph.local_supremal) > 1)
    {
        for (size_t i = 0; i < array_length(graph.local_supremal) - 1; ++i)
        {
            GVA_Variant variant;
            gva_edges(graph.observed.str,
                graph.local_supremal[i], graph.local_supremal[i + 1],
                i == 0, i == array_length(graph.local_supremal) - 2,
                &variant);
            fprintf(stream, "        \"" GVA_VARIANT_FMT "\"%s\n", GVA_VARIANT_PRINT(variant), i < array_length(graph.local_supremal) - 2 ? "," : "");
        } // for
    } // if
    fprintf(stream, "    ],\n    \"supremal\": \"" GVA_VARIANT_FMT "\",\n", GVA_VARIANT_PRINT(graph.supremal));
    fprintf(stream, "    \"canonical\": [\n");
    GVA_Variant* canonical = gva_canonical(gva_std_allocator, graph);
    for (size_t i = 0; i < array_length(canonical); ++i)
    {
        fprintf(stream, "        \"" GVA_VARIANT_FMT "\"%s\n", GVA_VARIANT_PRINT(canonical[i]), i < array_length(canonical) - 1 ? "," : "");
    } // for
    canonical = ARRAY_DESTROY(gva_std_allocator, canonical);

    fprintf(stream, "    ],\n    \"distance\": %u\n", graph.distance);
    fprintf(stream, "}\n");
} // lcs_graph_json


static void
lcs_graph_raw(FILE* const stream, GVA_LCS_Graph const graph)
{
    fprintf(stream, "#nodes: %zu\n#edges: %zu\n", array_length(graph.nodes), array_length(graph.edges));
    fprintf(stream, "distance: %u\n", graph.distance);
    fprintf(stream, "source: %u\n", graph.source);
    fprintf(stream, "local supremal:\n");
    if (array_length(graph.local_supremal) > 1)
    {
        for (size_t i = 0; i < array_length(graph.local_supremal) - 1; ++i)
        {
            GVA_Variant variant;
            gva_edges(graph.observed.str,
                graph.local_supremal[i], graph.local_supremal[i + 1],
                i == 0, i == array_length(graph.local_supremal) - 2,
                &variant);
            fprintf(stream, "    (%u, %u, %u) -> (%u, %u, %u): " GVA_VARIANT_FMT "\n",
                graph.local_supremal[i].row, graph.local_supremal[i].col, graph.local_supremal[i].length,
                graph.local_supremal[i + 1].row, graph.local_supremal[i + 1].col, graph.local_supremal[i + 1].length,
                GVA_VARIANT_PRINT(variant));
        } // for
    } // if

    if (graph.nodes == NULL)
    {
        return;
    } // if

    for (size_t i = 0; i < array_length(graph.nodes); ++i)
    {
        fprintf(stream, "%zu: (%u, %u, %u):\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length);
        if (graph.nodes[i].lambda != GVA_NULL)
        {
            fprintf(stream, "    lambda: (%u, %u, %u)\n", graph.nodes[graph.nodes[i].lambda].row, graph.nodes[graph.nodes[i].lambda].col, graph.nodes[graph.nodes[i].lambda].length);
        } // if
        for (gva_uint j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
        {
            fprintf(stream, "    %u: (%u, %u, %u): ", graph.edges[j].tail, graph.nodes[graph.edges[j].tail].row, graph.nodes[graph.edges[j].tail].col, graph.nodes[graph.edges[j].tail].length);
            GVA_Variant variant;
            gva_uint const count = gva_edges(graph.observed.str,
                graph.nodes[i], graph.nodes[graph.edges[j].tail],
                i == graph.source, graph.nodes[graph.edges[j].tail].edges == GVA_NULL,
                &variant);
            fprintf(stream, GVA_VARIANT_FMT " x %u\n", GVA_VARIANT_PRINT(variant), count);
        } // for
    } // for
} // lcs_graph_raw


int
extract(int const argc, char* argv[static argc + 1])
{
    if (argc != 3)
    {
        fprintf(stderr, "usage %s reference observed\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    /*
    char const* const restrict reference = "ATTCTATCTTCTGTCTACATAAGATGTCATACTAGAGGGCATATCTGCAATGTATACATATTATCTTTTCCAGCATGCATTCAGTTGTGTTGGAATAATTTATGTACACCTTTATAAACGCTGAGCCTCACAAGAGCCATGTGCCACGTATTGTTTTCTTACTACTTTTTGGGATACCTGGCACGTAATAGACACTCATTGAAAGTTTCCTAATGAATGAAGTACAAAGATAAAACAAGTTATAGACTGATTCTTTTGAGCTGTCAAGGTTGTAAATAGACTTTTGCTCAATCAATTCAAATGGTGGCAGGTAGTGGGGGTAGAGGGATTGGTATGAAAAACATAAGCTTTCAGAACTCCTGTGTTTATTTTTAGAATGTCAACTGCTTGAGTGTTTTTAACTCTGTGGTATCTGAACTATCTTCTCTAACTGCAGGTTGGGCTCAGATCTGTGATAGAACAGTTTCCTGGGAAGCTTGACTTTGTCCTTGTGGATGGGGGCTGTGTCCTAAGCCATGGCCACAAGCAGTTGATGTGCTTGGCTAGATCTGTTCTCAGTAAGGCGAAGATCTTGCTGCTTGATGAACCCAGTGCTCATTTGGATCCAGTGTGAGTTTCAGATGTTCTGTTACTTAATAGCACAGTGGGAACAGAATCATTATGCCTGCTTCATGGTGACACATATTTCTATTAGGCTGTCATGTCTGCGTGTGGGGGTCTCCCCCAAGATATGAAATAATTGCCCAGTGGAAATGAGCATAAATGCATATTTCCTTGCTAAGAGTCTTGTGTTTTCTTCCGAAGATAGTTTTTAGTTTCATACAAACTCTTCCCCCTTGTCAACACATGATGAAGCTTTTAAATACATGGGCCTAATCTGATCCTTATGATTTGCCTTTGTATCCCATTTATACCATAAGCATGTTTATAGCCCCAAATAAAGAAGTACTGGTGATTCTACATAATGAAAAATGTACTCATTTATTAAAGTTTCTTTGAAATATTTGTCCTGTTTATTTATGGATACTTAGAGTCTACCCCATGGTTGAAAAGCTGATTGTGGCTAACGCTATATCAACATTATGTGAAAAGAACTTAAAGAAATAAGTAATTTAAAGAGATAATAGAACAATAGACATATTATCAAGGTAAATACAGATCATTACTGTTCTGTGATATTATGTGTGGTATTTTCTTTCTTTTCTAGAACATACCAAATAATTAGAAGAACTCTAAAACAAGCATTTGCTGATTGCACAGTAATTCTCTGTGAACACAGGATAGAAGCAATGCTGGAATGCCAACAATTTTTGGTGAGTCTTTATAACTTTACTTAAGATCTCATTGCCCTTGTAATTCTTGATAACAATCTCACATGTGATAGTTCCTGCAAATTGCAACAATGTACAAGTTCTTTTCAAAAATATGTATCATACAGCCATCCAGCTTTACTCAAAATAGCTGCACAAGTTTTTCACTTTGATCTGAGCCATGTGGTGAGGTTGAAATATAGTAAATCTAAAATGGCAGCATATTACTAAGTTATGTTTATAAATAGGATATATATACTTTTTGAGCCCTTTATTTGGGGACCAAGTCATACAAAATACTCTACTGTTTAAGATTTTAAAAAAGGTCCCTGTGATTCTTTCAATAACTAAATGTCCCATGGATGTGGTCTGGGACAGGCCTAGTTGTCTTACAGTCTGATTTATGGTATTAATGACAAAGTTGAGAGGCACATTTCATTTTT";
    char const* const restrict observed = "ATTCTATCTTCTGTCTACATAAGATGTCATACTAGAGGGCATATCTGCAATGTATACATATTATCTTTTCCAGCATGCATTCAGTTGTGTTGGAATAATTTATGTACACCTTTATAAACGCTGAGCCTCACAAGAGCCATGTGCCACGTATTGTTTTCTTACTACTTTTTGGGATACCTGGCACGTAATAGACACTCATTGAAAGTTTCCTAATGAATGAAGTACAAAGATAAAACAAGTTATAGACTGATTCTTTTGAGCTGTCAAGGTTGTAAATAGACTTTTGCTCAATCAATTCAAATGGTGGCAGGTAGTGGGGGTAGAGGGATTGGTATGAAAAACATAAGCTTTCAGAACTCCTGTGTTTATTTTTAGAATGTCAACTGCTTGAGTGTTTTTAACTCTGTGGTATCTGAACTATCTTCTCTAACTGCAGGTGAGTCTTTATAACTTTACTTAAGATCTCATTGCCCTTGTAATTCTTGATAACAATCTCACATGTGATAGTTCCTGCAAATTGCAACAATGTACAAGTTCTTTTCAAAAATATGTATCATACAGCCATCCAGCTTTACTCAAAATAGCTGCACAAGTTTTTCACTTTGATCTGAGCCATGTGGTGAGGTTGAAATATAGTAAATCTAAAATGGCAGCATATTACTAAGTTATGTTTATAAATAGGATATATATACTTTTTGAGCCCTTTATTTGGGGACCAAGTCATACAAAATACTCTACTGTTTAAGATTTTAAAAAAGGTCCCTGTGATTCTTTCAATAACTAAATGTCCCATGGATGTGGTCTGGGACAGGCCTAGTTGTCTTACAGTCTGATTTATGGTATTAATGACAAAGTTGAGAGGCACATTTCATTTTT";
    */
    char const* const restrict reference = argv[1];
    char const* const restrict observed = argv[2];

    size_t const len_ref = strlen(reference);
    size_t const len_obs = strlen(observed);

    GVA_LCS_Graph graph = gva_lcs_graph_init(gva_std_allocator, len_ref, reference, len_obs, observed, 0);

    lcs_graph_raw(stderr, graph);
    lcs_graph_dot(stderr, graph);
    lcs_graph_json(stdout, graph);
    lcs_graph_svg(stderr, graph);

    GVA_Variant* canonical = gva_canonical(gva_std_allocator, graph);
    for (size_t i = 0; i < array_length(canonical); ++i)
    {
        fprintf(stderr, GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(canonical[i]));
    } // for
    canonical = ARRAY_DESTROY(gva_std_allocator, canonical);

    gva_lcs_graph_destroy(gva_std_allocator, graph);

    return EXIT_SUCCESS;
} // extract


void
index_dot(FILE* const stream, GVA_Stabbing_Index const index)
{
    fprintf(stream, "strict digraph{\nrankdir=BT\n");
    for (gva_uint i = 0; i < array_length(index.entries); ++i)
    {
        fprintf(stream, "%u[label=\"%u\\n[%u, %u]\"]\n", i, i, index.entries[i].start, index.entries[i].end);
        if (index.entries[i].parent != GVA_NULL)
        {
            fprintf(stream, "%u->%u\n", i, index.entries[i].parent);
        } // if
    } // for
    fprintf(stream, "}\n");
} // index_dot


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


void
trie_dot(FILE* const stream, Trie const self)
{
    fprintf(stream, "strict digraph{\n\"root\"[label=\"\",shape=point]\n");
    if (self.nodes != NULL)
    {
        for (gva_uint i = self.root; i != GVA_NULL; i = self.nodes[i].next)
        {
            int const p_len = self.nodes[i].end - self.nodes[i].p_start;
            fprintf(stream, "root->%u[label=\"%.*s\"]\n", i, p_len, self.strings.str + self.nodes[i].p_start);
        } // for
        for (size_t i = 0; i < array_length(self.nodes); ++i)
        {
            int const len = self.nodes[i].end - self.nodes[i].start;
            fprintf(stream, "%zu[label=\"%zu\\n%.*s\"]\n", i, i, len, self.strings.str + self.nodes[i].start);
            for (gva_uint j = self.nodes[i].link; j != GVA_NULL; j = self.nodes[j].next)
            {
                int const p_len = self.nodes[j].end - self.nodes[j].p_start;
                fprintf(stream, "%zu->%u[label=\"%.*s\"]\n", i, j, p_len, self.strings.str + self.nodes[j].p_start);
            } // for
        } // for
    } // if
    fprintf(stream, "}\n");
} // trie_dot


void
trie_raw(FILE* const stream, Trie const self)
{
    fprintf(stream, "root: %u\n", self.root);
    fprintf(stream, "strings (%zu): %.*s\n", self.strings.len, (int) self.strings.len, self.strings.str);
    fprintf(stream, "nodes (%zu):\n", array_length(self.nodes));
    for (size_t i = 0; i < array_length(self.nodes); ++i)
    {
        fprintf(stream, "[%zu]\n    .link: %2d\n    .next: %2d\n    .p_start: %2u\n    .start: %2u\n    .end:   %2u\n", i, self.nodes[i].link, self.nodes[i].next, self.nodes[i].p_start, self.nodes[i].start, self.nodes[i].end);
    } // for
} // trie_raw


void
bitset_print(size_t const bitset[static 1], size_t const size)
{
    for (intmax_t i = array_length(bitset) * sizeof(*bitset) * CHAR_BIT - 1; i >= 0; --i)
    {
        if ((size_t) i > size)
        {
            fprintf(stderr, " ");
        } // if
        else
        {
            fprintf(stderr, "%zu", bitset[i / (sizeof(*bitset) * CHAR_BIT)] >> (i % (sizeof(*bitset) * CHAR_BIT)) & 0x1);
        } // else
        if (i % (sizeof(*bitset) * CHAR_BIT) == 0)
        {
            fprintf(stderr, "\n");
        } // if
        else if (i % CHAR_BIT == 0)
        {
            fprintf(stderr, " ");
        } // if
    } // for
} // bitset_print


static inline int
entry_cmp(void const* lhs, void const* rhs)
{
    GVA_Stabbing_Entry const* const lhs_entry = lhs;
    GVA_Stabbing_Entry const* const rhs_entry = rhs;

    if (lhs_entry->start < rhs_entry->start)
    {
        return -1;
    } // if
    if (lhs_entry->start > rhs_entry->start)
    {
        return 1;
    } // if
    if (lhs_entry->end > rhs_entry->end)
    {
        return -1;
    } // if
    if (lhs_entry->end < rhs_entry->end)
    {
        return 1;
    } // if
    if (lhs_entry->distance > rhs_entry->distance)
    {
        return -1;
    } // if
    if (lhs_entry->distance < rhs_entry->distance)
    {
        return 1;
    } // if
    return 0;
} // entry_cmp


static size_t
binary_search(size_t const n, GVA_Stabbing_Entry const entries[static n],
    GVA_Stabbing_Entry const key, size_t low)
{
    size_t high = n;
    while (low <= high)
    {
        size_t const mid = low + (high - low) / 2;
        int const cmp = entry_cmp(&entries[mid], &key);
        if (cmp < 0)
        {
            low = mid + 1;
        } // if
        else if (cmp > 0 && mid > 0)
        {
            high = mid - 1;
        } // if
        else
        {
            return mid;
        } // else
    } // while
    return low;
} // binary_search


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
}


int
all(int argc, char* argv[static argc + 1])
{
    if (argc < 2)
    {
        fprintf(stderr, "usage: %s path_to_reference_fasta\n", argv[0]);
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
    reference = gva_fasta_sequence(gva_std_allocator, stream);
    fclose(stream);

    fprintf(stderr, "reference length: %zu\n", reference.len);

    TIC(perf_calculate_supremals);

    GVA_Stabbing_Index index = gva_stabbing_index_init(gva_std_allocator, reference.len);
    Trie trie = {{0, NULL}, NULL, GVA_NULL};

    size_t count = 0;
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
            fprintf(stderr, "error: SPDI parsing failed at line %zu: %s", count + 1, line);
            continue;
        } // if

        GVA_LCS_Graph graph = gva_lcs_graph_from_variants(gva_std_allocator, reference.len, reference.str, 1, &variant);

        size_t const allele = gva_stabbing_index_add_allele(gva_std_allocator, &index, rsid, graph.distance) - 1;
        for (size_t i = 0; i < array_length(graph.local_supremal) - 1; ++i)
        {
            GVA_Variant variant;
            gva_edges(graph.observed.str,
                graph.local_supremal[i], graph.local_supremal[i + 1],
                i == 0, i == array_length(graph.local_supremal) - 2,
                &variant);

            gva_uint const inserted = trie_insert(gva_std_allocator, &trie, variant.sequence.len, variant.sequence.str);
            gva_stabbing_index_add_part(gva_std_allocator, &index,
                variant.start, variant.end, inserted, graph.local_supremal[i + 1].edges, allele);

        } // for

        gva_string_destroy(gva_std_allocator, graph.observed);
        gva_lcs_graph_destroy(gva_std_allocator, graph);

        count += 1;
    } // while
    fprintf(stderr, "number of variants: %zu\n", array_length(index.entries) - 1);
    fprintf(stderr, "number of alleles: %zu\n", array_length(index.alleles));

    TOC(perf_calculate_supremals);

    TIC(perf_build_index);

    qsort(index.entries, array_length(index.entries), sizeof(*index.entries), entry_cmp);
    gva_stabbing_index_build(gva_std_allocator, index);

    for (size_t i = array_length(index.entries) - 1; i > 0; --i)
    {
        index.entries[i].next = index.alleles[index.entries[i].allele].head;
        index.alleles[index.entries[i].allele].head = i;
    } // for

    TOC(perf_build_index);

    for (size_t i = 0; i < array_length(index.alleles); ++i)
    {
        // fprintf(stderr, "%zu\n", index.alleles[i].data);
        for (gva_uint j = index.alleles[i].head; j != GVA_NULL; j = index.entries[j].next)
        {
            GVA_Variant const variant = {index.entries[j].start, index.entries[j].end, trie_string(trie, index.entries[j].inserted)};
            fprintf(stderr, GVA_VARIANT_FMT " ", GVA_VARIANT_PRINT(variant));
        } // for
        fprintf(stderr, "\n");
    } // for

    TIC(perf_calculate_relations);

    // NG_008376.4:9199:5:CCCCG
    // GVA_Variant rhs = {9199, 9203, {4, "CCCC"}};
    GVA_Variant rhs = {9199, 9204, {5, "CCCCG"}};
    gva_uint rhs_dist = 4;
    gva_uint* results = gva_stabbing_index_intersect(gva_std_allocator, index, rhs.start, rhs.end);

    for (size_t i = 0; i < array_length(results); ++i)
    {
        if (results[i] != 3858)
        {
            //continue;
        }
        GVA_Variant const lhs = {index.entries[results[i]].start, index.entries[results[i]].end, trie_string(trie, index.entries[results[i]].inserted)};
        GVA_Relation relation = compare_from_index(reference, lhs, index.entries[results[i]].distance, rhs, rhs_dist);

        // fprintf(stderr, GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(lhs));
        // fprintf(stderr, GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(rhs));

        // fprintf(stderr, "%u %u\n", index.entries[results[i]].distance, rhs_dist);
        if (relation == GVA_EQUIVALENT)
        {
            if (index.alleles[index.entries[results[i]].allele].distance > rhs_dist)
            {
                relation = GVA_CONTAINS;
            }
        }
        else if (relation == GVA_IS_CONTAINED)
        {
            if (index.alleles[index.entries[results[i]].allele].distance > index.entries[results[i]].distance)
            {
                relation = GVA_OVERLAP;
            }
        }

        printf("%u %zu %s\n", results[i], index.alleles[index.entries[results[i]].allele].data, GVA_RELATION_LABELS[relation]);
    }

    // printf("%s\n", GVA_RELATION_LABELS[relation]);

    results = ARRAY_DESTROY(gva_std_allocator, results);

    TOC(perf_calculate_relations);

    gva_stabbing_index_destroy(gva_std_allocator, &index);
    trie_destroy(gva_std_allocator, &trie);
    gva_string_destroy(gva_std_allocator, reference);
    return EXIT_SUCCESS;
} // all


int
main(int argc, char* argv[static argc + 1])
{
    (void) argv;
    //return extract(argc, argv);

    // all
    return all(argc, argv);
} // main
