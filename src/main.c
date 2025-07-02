#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // FILE, stderr, fclose, fopen, fprintf
#include <stdlib.h>     // EXIT_*
#include <string.h>     // strcmp, strlen


#include "../include/compare.h"     // gva_compare
#include "../include/edit.h"        // gva_edit_distance
#include "../include/extract.h"     // gva_extract
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, GVA_Variant, gva_lcs_graph_*, gva_edges
#include "../include/relations.h"   // GVA_RELATION_LABELS
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/types.h"       // GVA_NULL, GVA_String, gva_uint
#include "../include/utils.h"       // gva_fasta_sequence
#include "../include/variant.h"     // GVA_VARIANT_*, GVA_Variant, gva_parse_spdi
#include "array.h"      // ARRAY_DESTROY, array_length
#include "bitset.h"     // bitset_*


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
compare(void)
{
    GVA_String seq = {0, NULL};
    FILE* const restrict stream = fopen("data/NC_000001.11.fasta", "r");
    if (stream != NULL)
    {
        seq = gva_fasta_sequence(gva_std_allocator, stream);
        fclose(stream);
    } // if

    // NC_000001.11:169556314:6:TTTTTCTCT 3 NC_000001.11:169556314:8:TTTTTGTCTGGCTGC 7 4 is_contained
    char lhs_spdi[4096];
    size_t lhs_dist;
    char rhs_spdi[4096];
    size_t rhs_dist;
    size_t dist;
    char python[32];
    size_t count = 0;
    while (fscanf(stdin, "%4096s %zu %4096s %zu %zu %32s\n", lhs_spdi, &lhs_dist, rhs_spdi, &rhs_dist, &dist, python) == 6)
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
        count += 1;
//        printf("%d\n", gva_compare(gva_std_allocator, seq.len, seq.str, lhs, rhs));

        if (strcmp(python, GVA_RELATION_LABELS[gva_compare(gva_std_allocator, seq.len, seq.str, lhs, rhs)]) != 0)
        {
            printf("different %s %zu %s %zu %zu %s %s\n", lhs_spdi, lhs_dist, rhs_spdi, rhs_dist, dist, python, GVA_RELATION_LABELS[gva_compare(gva_std_allocator, seq.len, seq.str, lhs, rhs)]);
        }
//        printf("%s %s\n", lhs_spdi, rhs_spdi);

    }
    printf("%zu\n", count);
    seq.str = gva_std_allocator.allocate(gva_std_allocator.context, seq.str, seq.len, 0);

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
    size_t rsid;
    char spdi[4096];
    char hgvs[4096];
    char python_sup[4096];
    //char c_sup[4096];
    size_t python_distance;
    size_t python_nodes;
    size_t python_edges;
    size_t count = 0;
    size_t counter = 0;
    while (fscanf(stdin, "%zu %4096s %4096s %4096s %zu %zu %zu\n", &rsid, spdi, hgvs, python_sup, &python_distance, &python_nodes, &python_edges) == 7)
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

        gva_uint sup_len = graph.supremal.end - graph.supremal.start;
        size_t* dels = bitset_init(gva_std_allocator, sup_len + 1);
        size_t* as = bitset_init(gva_std_allocator, sup_len + 1);
        size_t* cs = bitset_init(gva_std_allocator, sup_len + 1);
        size_t* gs = bitset_init(gva_std_allocator, sup_len + 1);
        size_t* ts = bitset_init(gva_std_allocator, sup_len + 1);

        counter += bitset_fill(graph, graph.supremal.start, graph.supremal.start, graph.supremal.end, dels, as, cs, gs, ts);

        //for (size_t i = 0; i < array_length(dels); ++i)
        //{
        //    fprintf(stderr, "%2zu: %016zx\n", i, dels[i]);
        //    fprintf(stderr, "%2zu: %016zx\n", i, as[i]);
        //    fprintf(stderr, "%2zu: %016zx\n", i, cs[i]);
        //    fprintf(stderr, "%2zu: %016zx\n", i, gs[i]);
        //    fprintf(stderr, "%2zu: %016zx\n", i, ts[i]);
        //} // for

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


int
main(int argc, char* argv[static argc + 1])
{
    (void) argv;
    // check python relation output
    //return compare();

    // calculate graphs for all inputs
    return all_graphs();

    //extract single variant
    //return extract(argc, argv);
} // main
