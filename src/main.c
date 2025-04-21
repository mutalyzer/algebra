#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint32_t
#include <stdio.h>      // stderr, fprintf, printf
#include <stdlib.h>     // EXIT_*, atoi, rand, srand
#include <string.h>     // strlen
#include "../include/array.h"           // va_array_*
#include "../include/edit.h"            // va_edit
#include "../include/graph.h"           // Graph, to_dot
#include "../include/graph2.h"          // build
#include "../include/std_alloc.h"       // va_std_allocator


#define print_variant(variant, observed) variant.start, variant.end, (int) variant.obs_end - variant.obs_start, observed + variant.obs_start
#define VAR_FMT "%u:%u/%.*s"


static uint32_t const GVA_NULL = UINT32_MAX;


static size_t
random_sequence(size_t const min, size_t const max, char sequence[static max])
{
    size_t const len = rand() % (max - min) + min;

    for (size_t i = 0; i < len; ++i)
    {
        sequence[i] = "ACGT"[rand() % 4];
    } // for

    return len;
} // random_sequence


static int
check(size_t const len_ref, char const reference[static len_ref],
      size_t const len_obs, char const observed[static len_obs],
      bool const debug)
{
    VA_LCS_Node** lcs_nodes = NULL;
    size_t const len_lcs = va_edit(va_std_allocator, len_ref, reference, len_obs, observed, &lcs_nodes);

    if (debug)
    {
        fprintf(stderr, "%zu\n", len_lcs);
        for (size_t i = 0; i < len_lcs; ++i)
        {
            fprintf(stderr, "%zu:  ", i);
            for (size_t j = 0; j < va_array_length(lcs_nodes[i]); ++j)
            {
                fprintf(stderr, "(%u, %u, %u) ", lcs_nodes[i][j].row, lcs_nodes[i][j].col, lcs_nodes[i][j].length);
            } // for
            fprintf(stderr, "\n");
        } // for
    } // if

    Graph graph = build_graph(va_std_allocator, len_ref, len_obs, len_lcs, lcs_nodes, 0, false);

    if (debug)
    {
        to_dot(graph, len_obs, observed);
    } // if

    printf("#nodes: %zu\n#edges: %zu\n", va_array_length(graph.nodes), va_array_length(graph.edges));
    printf("source: %u\n", graph.source);
    for (size_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        printf("%zu: (%u, %u, %u):\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length);
        if (graph.nodes[i].lambda != GVA_NULL)
        {
            printf("    (%u, %u, %u): lambda\n", graph.nodes[graph.nodes[i].lambda].row, graph.nodes[graph.nodes[i].lambda].col, graph.nodes[graph.nodes[i].lambda].length);
        }
        for (uint32_t j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
        {
            printf("    (%u, %u, %u): " VAR_FMT "\n", graph.nodes[graph.edges[j].tail].row, graph.nodes[graph.edges[j].tail].col, graph.nodes[graph.edges[j].tail].length, print_variant(graph.edges[j].variant, observed));
        } // for
    } // for

    size_t count = 0;
    for (uint32_t head = 0; head < va_array_length(graph.nodes); ++head)
    {
        uint32_t* tails = va_array_init(va_std_allocator, 10, sizeof(*tails));
        for (uint32_t i = graph.nodes[head].edges; i != (uint32_t) -1; i = graph.edges[i].next)
        {
            uint32_t const tail = graph.edges[i].tail;

            bool found = false;
            for (uint32_t j = 0; j < va_array_length(tails); ++j)
            {
                if (tail == tails[j])
                {
                    found = true;
                    break;
                } // if
            } // for

            va_array_append(va_std_allocator, tails, tail);

            size_t const n = edges(
                graph.nodes[head].row, graph.nodes[head].col, graph.nodes[head].length,
                graph.nodes[tail].row, graph.nodes[tail].col, graph.nodes[tail].length + (graph.nodes[tail].edges == (uint32_t) -1),
                graph.source == head,
                len_obs, observed,
                graph.edges[i].variant
            );
            if (n == (size_t) -1)
            {
                printf("EDGE NOT FOUND: " VAR_FMT "\n", print_variant(graph.edges[i].variant, observed));
                return EXIT_FAILURE;
            } // if
            if (!found)
            {
                count += n;
            } // if
        } // for
        tails = va_array_destroy(va_std_allocator, tails);
    } // for

    if (count != va_array_length(graph.edges))
    {
        printf("COUNTS: %zu vs %zu\n", count, va_array_length(graph.edges));
        return EXIT_FAILURE;
    } // if

    destroy(va_std_allocator, &graph);
    return EXIT_SUCCESS;
} // check


static void
print_graph(Graph const graph, size_t const len_obs, char const observed[static len_obs]) {
    fprintf(stderr, "#nodes: %zu\n#edges: %zu\n", va_array_length(graph.nodes), va_array_length(graph.edges));
    fprintf(stderr, "source: %u\n", graph.source);
    if (graph.nodes == NULL)
    {
        return;
    } // if
    for (size_t i = 0; i < va_array_length(graph.nodes); ++i) {
        fprintf(stderr, "%zu: (%u, %u, %u):\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length);
        if (graph.nodes[i].lambda != GVA_NULL) {
            fprintf(stderr, "    (%u, %u, %u): lambda\n", graph.nodes[graph.nodes[i].lambda].row,
                    graph.nodes[graph.nodes[i].lambda].col, graph.nodes[graph.nodes[i].lambda].length);
        } // if
        for (uint32_t j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next) {
            fprintf(stderr, "    (%u, %u, %u): ", graph.nodes[graph.edges[j].tail].row,
                    graph.nodes[graph.edges[j].tail].col, graph.nodes[graph.edges[j].tail].length);
            fprintf(stderr, VAR_FMT "\n", print_variant(graph.edges[j].variant, observed));
        } // for
    } // for
} // print_graph


int
main(int argc, char* argv[static argc + 1])
{
/*
    (void) argv;
    char const* const restrict reference = "ATTCTATCTTCTGTCTACATAAGATGTCATACTAGAGGGCATATCTGCAATGTATACATATTATCTTTTCCAGCATGCATTCAGTTGTGTTGGAATAATTTATGTACACCTTTATAAACGCTGAGCCTCACAAGAGCCATGTGCCACGTATTGTTTTCTTACTACTTTTTGGGATACCTGGCACGTAATAGACACTCATTGAAAGTTTCCTAATGAATGAAGTACAAAGATAAAACAAGTTATAGACTGATTCTTTTGAGCTGTCAAGGTTGTAAATAGACTTTTGCTCAATCAATTCAAATGGTGGCAGGTAGTGGGGGTAGAGGGATTGGTATGAAAAACATAAGCTTTCAGAACTCCTGTGTTTATTTTTAGAATGTCAACTGCTTGAGTGTTTTTAACTCTGTGGTATCTGAACTATCTTCTCTAACTGCAGGTTGGGCTCAGATCTGTGATAGAACAGTTTCCTGGGAAGCTTGACTTTGTCCTTGTGGATGGGGGCTGTGTCCTAAGCCATGGCCACAAGCAGTTGATGTGCTTGGCTAGATCTGTTCTCAGTAAGGCGAAGATCTTGCTGCTTGATGAACCCAGTGCTCATTTGGATCCAGTGTGAGTTTCAGATGTTCTGTTACTTAATAGCACAGTGGGAACAGAATCATTATGCCTGCTTCATGGTGACACATATTTCTATTAGGCTGTCATGTCTGCGTGTGGGGGTCTCCCCCAAGATATGAAATAATTGCCCAGTGGAAATGAGCATAAATGCATATTTCCTTGCTAAGAGTCTTGTGTTTTCTTCCGAAGATAGTTTTTAGTTTCATACAAACTCTTCCCCCTTGTCAACACATGATGAAGCTTTTAAATACATGGGCCTAATCTGATCCTTATGATTTGCCTTTGTATCCCATTTATACCATAAGCATGTTTATAGCCCCAAATAAAGAAGTACTGGTGATTCTACATAATGAAAAATGTACTCATTTATTAAAGTTTCTTTGAAATATTTGTCCTGTTTATTTATGGATACTTAGAGTCTACCCCATGGTTGAAAAGCTGATTGTGGCTAACGCTATATCAACATTATGTGAAAAGAACTTAAAGAAATAAGTAATTTAAAGAGATAATAGAACAATAGACATATTATCAAGGTAAATACAGATCATTACTGTTCTGTGATATTATGTGTGGTATTTTCTTTCTTTTCTAGAACATACCAAATAATTAGAAGAACTCTAAAACAAGCATTTGCTGATTGCACAGTAATTCTCTGTGAACACAGGATAGAAGCAATGCTGGAATGCCAACAATTTTTGGTGAGTCTTTATAACTTTACTTAAGATCTCATTGCCCTTGTAATTCTTGATAACAATCTCACATGTGATAGTTCCTGCAAATTGCAACAATGTACAAGTTCTTTTCAAAAATATGTATCATACAGCCATCCAGCTTTACTCAAAATAGCTGCACAAGTTTTTCACTTTGATCTGAGCCATGTGGTGAGGTTGAAATATAGTAAATCTAAAATGGCAGCATATTACTAAGTTATGTTTATAAATAGGATATATATACTTTTTGAGCCCTTTATTTGGGGACCAAGTCATACAAAATACTCTACTGTTTAAGATTTTAAAAAAGGTCCCTGTGATTCTTTCAATAACTAAATGTCCCATGGATGTGGTCTGGGACAGGCCTAGTTGTCTTACAGTCTGATTTATGGTATTAATGACAAAGTTGAGAGGCACATTTCATTTTT";
    char const* const restrict observed = "ATTCTATCTTCTGTCTACATAAGATGTCATACTAGAGGGCATATCTGCAATGTATACATATTATCTTTTCCAGCATGCATTCAGTTGTGTTGGAATAATTTATGTACACCTTTATAAACGCTGAGCCTCACAAGAGCCATGTGCCACGTATTGTTTTCTTACTACTTTTTGGGATACCTGGCACGTAATAGACACTCATTGAAAGTTTCCTAATGAATGAAGTACAAAGATAAAACAAGTTATAGACTGATTCTTTTGAGCTGTCAAGGTTGTAAATAGACTTTTGCTCAATCAATTCAAATGGTGGCAGGTAGTGGGGGTAGAGGGATTGGTATGAAAAACATAAGCTTTCAGAACTCCTGTGTTTATTTTTAGAATGTCAACTGCTTGAGTGTTTTTAACTCTGTGGTATCTGAACTATCTTCTCTAACTGCAGGTGAGTCTTTATAACTTTACTTAAGATCTCATTGCCCTTGTAATTCTTGATAACAATCTCACATGTGATAGTTCCTGCAAATTGCAACAATGTACAAGTTCTTTTCAAAAATATGTATCATACAGCCATCCAGCTTTACTCAAAATAGCTGCACAAGTTTTTCACTTTGATCTGAGCCATGTGGTGAGGTTGAAATATAGTAAATCTAAAATGGCAGCATATTACTAAGTTATGTTTATAAATAGGATATATATACTTTTTGAGCCCTTTATTTGGGGACCAAGTCATACAAAATACTCTACTGTTTAAGATTTTAAAAAAGGTCCCTGTGATTCTTTCAATAACTAAATGTCCCATGGATGTGGTCTGGGACAGGCCTAGTTGTCTTACAGTCTGATTTATGGTATTAATGACAAAGTTGAGAGGCACATTTCATTTTT";
*/

    if (argc == 2)
    {
        srand(atoi(argv[1]));
        size_t count = 0;
        while (true)
        {
            char reference[4096];
            char observed[4096];
            size_t const len_ref = random_sequence(0, 40, reference);
            size_t const len_obs = random_sequence(0, 40, observed);

            printf("%zu: %.*s %.*s\n", count, (int) len_ref, reference, (int) len_obs, observed);
            if (check(len_ref, reference, len_obs, observed, false) == EXIT_FAILURE)
            {
                return EXIT_FAILURE;
            } // if
            count += 1;
        } // while
    } // if
    else if (argc == 3)
    {
        char const* const restrict reference = argv[1];
        char const* const restrict observed = argv[2];
        size_t const len_ref = strlen(reference);
        size_t const len_obs = strlen(observed);
        //check(len_ref, reference, len_obs, observed, true);
        Graph2 graph2 = build(len_ref, reference, len_obs, observed, 0);

        //to_json2(graph2, len_obs, observed);

        VA_LCS_Node** lcs_nodes = NULL;
        size_t const len_lcs = va_edit(va_std_allocator, len_ref, reference, len_obs, observed, &lcs_nodes);

        Graph graph = build_graph(va_std_allocator, len_ref, len_obs, len_lcs, lcs_nodes, 0, false);

        print_graph(graph, len_obs, observed);
        // to_dot(graph, len_obs, observed);
        // to_json(graph, len_obs, observed, false);

        // destroy graph2
        va_array_destroy(va_std_allocator, graph2.nodes);
        va_array_destroy(va_std_allocator, graph2.edges);

        destroy(va_std_allocator, &graph);

    } // if
    else
    {
        fprintf(stderr, "usage: %s\n", argv[0]);
        return EXIT_FAILURE;
    } // else

/*
    printf("nodes (%zu)\n  #\t(row, col, len)\tedges\tlambda\n", va_array_length(graph.nodes));
    for (size_t i = 0; i < va_array_length(graph.nodes); ++i)

    {
        printf("%c%2zu:\t(%u, %u, %u)\t%5d\t%6d\n", i == graph.source ? '*' : ' ', i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, (signed) graph.nodes[i].edges, (signed) graph.nodes[i].lambda);
    } // for
    printf("edges (%zu)\n  #\ttail\t\"start:end/seq\"\tnext\n", va_array_length(graph.edges));
    for (size_t i = 0; i < va_array_length(graph.edges); ++i)
    {
        printf("%3zu:\t%4u\t  \"%u:%u/%.*s\"\t%4d\n", i, graph.edges[i].tail, graph.edges[i].variant.start, graph.edges[i].variant.end, (int) graph.edges[i].variant.obs_end - graph.edges[i].variant.obs_start, observed + graph.edges[i].variant.obs_start, (signed) graph.edges[i].next);
    } // for
*/

    //local_supremal(va_std_allocator, graph, len_obs, observed);

    //canonical(va_std_allocator, graph, len_obs, observed);

    //to_json(graph, len_obs, observed, lambda);

    return EXIT_SUCCESS;
} // main
