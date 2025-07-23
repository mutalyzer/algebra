//#define _XOPEN_SOURCE
//#include <assert.h>
#include <errno.h>      // errno
#include <limits.h>     // CHAR_BIT
#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // FILE, stderr, fclose, fopen, fprintf
#include <stdlib.h>     // EXIT_*
#include <string.h>     // strcmp, strerror, strlen
//#include <unistd.h>     // optarg


#include "../include/compare.h"     // gva_compare
#include "../include/edit.h"        // gva_edit_distance
#include "../include/extract.h"     // gva_extract
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, GVA_Variant, gva_lcs_graph_*, gva_edges
#include "../include/relations.h"   // GVA_RELATION_LABELS
#include "../include/stabbing.h"    // GVA_Stabbing_*, gva_stabbing_*
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/types.h"       // GVA_NULL, GVA_String, gva_uint
#include "../include/utils.h"       // gva_fasta_sequence
#include "../include/variant.h"     // GVA_VARIANT_*, GVA_Variant, gva_parse_spdi
#include "array.h"      // ARRAY_DESTROY, array_length
#include "bitset.h"     // bitset_*
#include "trie.h"       // Trie, trie_*


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
    GVA_Variant* canonical = gva_extract(gva_std_allocator, graph);
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
compare(int argc, char* argv[static argc + 1])
{
    GVA_String reference = {0, NULL};
    errno = 0;
    FILE* const restrict stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if
    reference = gva_fasta_sequence(gva_std_allocator, stream);
    fclose(stream);

    disjoint_count = 0;

    // NC_000001.11:169556314:6:TTTTTCTCT 3 NC_000001.11:169556314:8:TTTTTGTCTGGCTGC 7 4 is_contained
    // NC_000001.10:g.103496806_103496807delAA NC_000001.10:103496805:13:AAAAAAAAAAA NC_000001.10:g.103496811A>C NC_000001.10:103496805:13:AAAAACAAAAAAA overlap
    char lhs_hgvs[4096] = {0};
    char lhs_spdi[4096] = {0};
    //size_t lhs_dist = 0;
    char rhs_hgvs[4096] = {0};
    char rhs_spdi[4096] = {0};
    //size_t rhs_dist = 0;
    //size_t dist = 0;
    //char python[32] = {0};
    char relation[32] = {0};
    size_t count = 0;
    //while (fscanf(stdin, "%4096s %zu %4096s %zu %zu %32s\n", lhs_spdi, &lhs_dist, rhs_spdi, &rhs_dist, &dist, python) == 6)
    while (fscanf(stdin, "%4095s %4095s %4095s %4095s %31s\n", lhs_hgvs, lhs_spdi, rhs_hgvs, rhs_spdi, relation) == 5)
    {
        GVA_Variant lhs;
        if (!gva_parse_spdi(strlen(lhs_spdi), lhs_spdi, &lhs))
        {
            fprintf(stderr, "PARSE ERROR: %s\n", lhs_spdi);
            continue;
        } // if
        GVA_Variant rhs;
        if (!gva_parse_spdi(strlen(rhs_spdi), rhs_spdi, &rhs))
        {
            fprintf(stderr, "PARSE ERROR: %s\n", rhs_spdi);
            continue;
        } // if

        //fprintf(stderr, "%s %zu %s %zu %zu %s\n", lhs_spdi, lhs_dist, rhs_spdi, rhs_dist, dist, python);
        size_t const relation = gva_compare(gva_std_allocator, reference.len, reference.str, lhs, rhs);
        //if (strcmp(python, GVA_RELATION_LABELS[relation]) != 0)
        {
            //printf("%s %zu %s %zu %zu %s %s\n", lhs_spdi, lhs_dist, rhs_spdi, rhs_dist, dist, python, GVA_RELATION_LABELS[relation]);
            printf("%s %s %s %s %s\n", lhs_hgvs, lhs_spdi, rhs_hgvs, rhs_spdi, GVA_RELATION_LABELS[relation]);
        } // if

        count += 1;
    } // while
    fprintf(stderr, "total: %zu\n", count);
    fprintf(stderr, "disjoint: %zu\n", disjoint_count);
    reference.str = gva_std_allocator.allocate(gva_std_allocator.context, reference.str, reference.len, 0);

    return 0;
} // compare


int
all_graphs(void)
{
    GVA_String seq = {0, NULL};
    FILE* const restrict stream = fopen("data/NC_000001.11.fasta", "r");
    if (stream != NULL)
    {
        seq = gva_fasta_sequence(gva_std_allocator, stream);
        fclose(stream);
    } // if

    //bool check = false;
    size_t rsid = 0;
    char spdi[4096] = {0};
    char hgvs[4096] = {0};
    char python_sup[4096] = {0};
    //char c_sup[4096] = {0};
    size_t python_distance = 0;
    size_t python_nodes = 0;
    size_t python_edges = 0;
    size_t count = 0;
    size_t counter = 0;
    while (fscanf(stdin, "%zu %4095s %4095s %4095s %zu %zu %zu\n", &rsid, spdi, hgvs, python_sup, &python_distance, &python_nodes, &python_edges) == 7)
    {
        GVA_Variant variant;
        if (!gva_parse_spdi(strlen(spdi), spdi, &variant))
        {
            fprintf(stderr, "PARSE ERROR: %zu %s %s %s %zu %zu %zu\n", rsid, spdi, hgvs, python_sup, python_distance, python_nodes, python_edges);
            continue;
        } // if
        count += 1;
        GVA_LCS_Graph graph = gva_lcs_graph_from_variant(gva_std_allocator, seq.len, seq.str, variant);

        //fprintf(stdout, "%zu %s %s %s %zu %zu %zu " GVA_VARIANT_FMT " %u %zu %zu\n",
        //        rsid, spdi, hgvs, python_sup, python_distance, python_nodes, python_edges, GVA_VARIANT_PRINT(graph.supremal), graph.distance, array_length(graph.nodes), array_length(graph.edges));
        //fprintf(stdout, "%zu %s %s %s %zu %zu %zu " GVA_VARIANT_FMT " %u %zu %zu\n",
        //        rsid, spdi, hgvs, python_sup, python_distance, python_nodes, python_edges, GVA_VARIANT_PRINT(graph.supremal), graph.distance, array_length(graph.nodes), array_length(graph.edges));
        //sprintf(c_sup, GVA_VARIANT_FMT, GVA_VARIANT_PRINT(graph.supremal));
        //fprintf(stderr, "%s %s\n", python_sup, c_sup);
        //if (check && strcmp(python_sup, c_sup))
        //{
        //    fprintf(stderr, "Different supremal!\n");
        //}

        //lcs_graph_raw(stderr, graph);
        //lcs_graph_dot(stderr, graph);

        //gva_uint const sup_len = graph.supremal.end - graph.supremal.start;
        //size_t* dels = bitset_init(gva_std_allocator, sup_len + 1);
        //size_t* as = bitset_init(gva_std_allocator, sup_len + 1);
        //size_t* cs = bitset_init(gva_std_allocator, sup_len + 1);
        //size_t* gs = bitset_init(gva_std_allocator, sup_len + 1);
        //size_t* ts = bitset_init(gva_std_allocator, sup_len + 1);

        //counter += bitset_fill(graph, graph.supremal.start, graph.supremal.start, graph.supremal.end, dels, as, cs, gs, ts);

        //ts = bitset_destroy(gva_std_allocator, ts);
        //gs = bitset_destroy(gva_std_allocator, gs);
        //cs = bitset_destroy(gva_std_allocator, cs);
        //as = bitset_destroy(gva_std_allocator, as);
        //dels = bitset_destroy(gva_std_allocator, dels);

        graph.observed.str = gva_std_allocator.allocate(gva_std_allocator.context, graph.observed.str, graph.observed.len, 0);
        gva_lcs_graph_destroy(gva_std_allocator, graph);
    } // while

    fprintf(stderr, "%zu %zu\n", count, counter);

    seq.str = gva_std_allocator.allocate(gva_std_allocator.context, seq.str, seq.len, 0);

    return EXIT_SUCCESS;
} // all_graphs


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

    GVA_Variant* canonical = gva_extract(gva_std_allocator, graph);
    for (size_t i = 0; i < array_length(canonical); ++i)
    {
        fprintf(stderr, GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(canonical[i]));
    } // for
    canonical = ARRAY_DESTROY(gva_std_allocator, canonical);

    /*
    size_t const distance = gva_edit_distance(gva_std_allocator, len_ref, reference, len_obs, observed);
    fprintf(stderr, "distance only: %zu\n", distance);
    */

    gva_lcs_graph_destroy(gva_std_allocator, graph);

    return EXIT_SUCCESS;
} // extract


static void
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


static inline bool
intersect(size_t const lhs_start, size_t const lhs_end,
    size_t const rhs_start, size_t const rhs_end)
{
    return rhs_end >= lhs_start && rhs_start <= lhs_end;
} // intersect


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


typedef struct
{
    gva_uint start;
    gva_uint end;
} Interval;


static inline int
interval_cmp(void const* lhs, void const* rhs)
{
    Interval const* const lhs_interval = lhs;
    Interval const* const rhs_interval = rhs;

    if (lhs_interval->start < rhs_interval->start)
    {
        return -1;
    } // if
    if (lhs_interval->start > rhs_interval->start)
    {
        return 1;
    } // if
    if (lhs_interval->end > rhs_interval->end)
    {
        return -1;
    } // if
    if (lhs_interval->end < rhs_interval->end)
    {
        return 1;
    } // if
    return 0;
} // interval cmp


int
faststabber(int const argc, char* argv[static argc + 1])
{
    if (argc < 2)
    {
        fprintf(stderr, "usage: %s path_to_variant_tsv\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    errno = 0;
    FILE* stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    //248956422
    //GVA_Stabbing_Index index = gva_stabbing_index_init(gva_std_allocator, 248956422);

    Interval* intervals = NULL;

    size_t line_count = 0;
    static char line[4096] = {0};
    while (fgets(line, sizeof(line), stream) != NULL)
    {
        size_t idx = 0;

        // ignore rsid
        parse_number(line, &idx);
        idx += 1;
        size_t const start = parse_number(line, &idx);
        idx += 1;
        size_t const end = parse_number(line, &idx);
        // ignore remainder of the line

        //gva_stabbing_index_add(gva_std_allocator, &index, start, end);
        ARRAY_APPEND(gva_std_allocator, intervals, ((Interval) {start, end}));
        line_count += 1;
    } // while

    errno = 0;
    if (feof(stream) == 0)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        fclose(stream);
        //gva_stabbing_index_destroy(gva_std_allocator, &index);
        intervals = ARRAY_DESTROY(gva_std_allocator, intervals);
        return EXIT_FAILURE;
    } // if
    fclose(stream);

    fprintf(stderr, "line count:  %zu\n", line_count);
    if (intervals == NULL)
    {
        return EXIT_FAILURE;
    } // if

    //qsort(intervals, array_length(intervals), sizeof(*intervals), interval_cmp);
    //fprintf(stderr, "entry count: %zu\n", array_length(index.entries) - 1);

    /*
    for (size_t i = 0; i < array_length(intervals); ++i)
    {
        printf("%u %u\n", intervals[i].start, intervals[i].end);
    } // for
    */

    //gva_stabbing_index_build(gva_std_allocator, index);

    //index_dot(stderr, index);
    //gva_uint* results = gva_stabbing_index_self_intersect(gva_std_allocator, index);
    //fprintf(stderr, "total: %zu\n", array_length(results));
    //results = ARRAY_DESTROY(gva_std_allocator, results);
/*
    size_t count = 0;
    for (size_t i = 1; i < array_length(index.entries); ++i)
    {
        gva_uint* results = gva_stabbing_index_intersect(gva_std_allocator, index, index.entries[i].start, index.entries[i].end);
        //count += array_length(results);
        for (size_t j = 0; j < array_length(results); ++j)
        {
            if (i < results[j])
            {
                //fprintf(stdout, "%2zu [%u, %u] %2u [%u, %u]\n", i, index.entries[i].start, index.entries[i].end, results[j], index.entries[results[j]].start, index.entries[results[j]].end);
                count += 1;
            } // if
        } // for
        results = ARRAY_DESTROY(gva_std_allocator, results);
    } // for
    fprintf(stderr, "total: %zu\n", count);
*/

    size_t count = 0;
    for (size_t i = 0; i < array_length(intervals); ++i)
    {
        for (size_t j = i + 1; j < array_length(intervals); ++j)
        {
            if (intervals[i].end < intervals[j].start)
            {
                break;
            } // if
            if (intersect(intervals[i].start, intervals[i].end, intervals[j].start, intervals[j].end))
            {
                //fprintf(stdout, "%2zu [%u, %u] %2zu [%u, %u]\n", i, index.entries[i].start, index.entries[i].end, j, index.entries[j].start, index.entries[j].end);
                count += 1;
            } // if
        } // for
    } // for
    fprintf(stderr, "total: %zu\n", count);

//    gva_stabbing_index_destroy(gva_std_allocator, &index);
    intervals = ARRAY_DESTROY(gva_std_allocator, intervals);

    return EXIT_SUCCESS;
} // faststabber


static void
trie_dot(FILE* const stream, Trie const self)
{
    fprintf(stream, "strict digraph{\n\"root\"[label=\"\",shape=point]\n");
    if (self.nodes != NULL)
    {
        for (gva_uint i = self.root; i != GVA_NULL; i = self.nodes[i].next)
        {
            int const p_len = self.nodes[i].end - self.nodes[i].p_start;
            fprintf(stream, "root->%u[label=\"%.*s\"]\n", i, p_len, self.strings + self.nodes[i].p_start);
        } // for
        for (size_t i = 0; i < array_length(self.nodes); ++i)
        {
            int const len = self.nodes[i].end - self.nodes[i].start;
            fprintf(stream, "%zu[label=\"%zu\\n%.*s\"]\n", i, i, len, self.strings + self.nodes[i].start);
            for (gva_uint j = self.nodes[i].link; j != GVA_NULL; j = self.nodes[j].next)
            {
                int const p_len = self.nodes[j].end - self.nodes[j].p_start;
                fprintf(stream, "%zu->%u[label=\"%.*s\"]\n", i, j, p_len, self.strings + self.nodes[j].p_start);
            } // for
        } // for
    } // if
    fprintf(stream, "}\n");
} // trie_dot


static void
trie_raw(FILE* const stream, Trie const self)
{
    fprintf(stream, "strings (%zu): %.*s\n", array_length(self.strings), (int) array_length(self.strings), self.strings);
    fprintf(stream, "nodes (%zu):\n", array_length(self.nodes));
    for (size_t i = 0; i < array_length(self.nodes); ++i)
    {
        fprintf(stream, "[%zu]\n    .link: %2d\n    .next: %2d\n    .p_start: %2u\n    .start: %2u\n    .end:   %2u\n", i, self.nodes[i].link, self.nodes[i].next, self.nodes[i].p_start, self.nodes[i].start, self.nodes[i].end);
    } // for
} // trie_raw


int
trie(int argc, char* argv[static argc + 1])
{
    if (argc < 2)
    {
        fprintf(stderr, "usage %s path_to_strings\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    errno = 0;
    FILE* stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return EXIT_FAILURE;
    } // if

    Trie trie = {NULL, NULL, GVA_NULL};

    size_t count = 0;
    static char line[4096] = {0};
    while (fgets(line, sizeof(line), stream) != NULL)
    {
        size_t const len = strlen(line) - 1;
        if (len == 1 && line[0] == '.')
        {
           continue;
        } // if
        line[len] = '\0';
        fprintf(stderr, "INSERT %.*s\n", (int) len, line);
        gva_uint const idx = trie_insert(gva_std_allocator, &trie, len, line);
        fprintf(stderr, "%u\n", idx);
        count += 1;
    } // while

    errno = 0;
    if (feof(stream) == 0)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        fclose(stream);
        trie_destroy(gva_std_allocator, &trie);
        return EXIT_FAILURE;
    } // if
    fclose(stream);

    fprintf(stderr, "total: %zu\n", count);

    //trie_print(stderr, trie, 0);
    //trie_raw(stderr, trie);
    trie_dot(stdout, trie);

    fprintf(stderr, "str len: %zu\n", array_length(trie.strings));
    fprintf(stderr, "#nodes : %zu\n", array_length(trie.nodes));
    //fprintf(stderr, "%p\n", (void*) trie_find(trie, 1, "AA"));


    trie.strings = ARRAY_DESTROY(gva_std_allocator, trie.strings);
    trie.nodes = ARRAY_DESTROY(gva_std_allocator, trie.nodes);

    trie_destroy(gva_std_allocator, &trie);

    return EXIT_SUCCESS;
} // trie


int
main(int argc, char* argv[static argc + 1])
{
    (void) argv;

    // check python relation output
    //return compare(argc, argv);

    // calculate graphs for all inputs
    //return all_graphs();

    // extract single variant
    //return extract(argc, argv);

    // faststabber
    //return faststabber(argc, argv);
    // trie
    return trie(argc, argv);
} // main
