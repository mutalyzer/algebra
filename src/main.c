#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint32_t
#include <stdio.h>      // stderr, fprintf, printf
#include <stdlib.h>     // EXIT_*
#include <string.h>     // strlen

#include "../include/alloc.h"           // VA_Allocator
#include "../include/array.h"           // va_array_*
#include "../include/bitset.h"          // VA_Bitset, va_bitset_*
#include "../include/edit.h"            // va_edit
#include "../include/queue.h"           // VA_Queue, va_queue_*
#include "../include/std_alloc.h"       // va_std_allocator
#include "../include/variant.h"         // VA_Variant


typedef struct
{
    uint32_t sink;
    uint32_t next;
    VA_Variant variant;
} Edge;


typedef struct
{
    uint32_t row;
    uint32_t col;
    uint32_t length;
    uint32_t edges;
    uint32_t lambda;
} Node;


typedef struct
{
    Node* nodes;
    Edge* edges;
    uint32_t source;
} Graph;


static inline size_t
add_node(VA_Allocator const allocator, Graph* const graph, size_t const row, size_t const col, size_t const length)
{
    va_array_append(allocator, graph->nodes, ((Node) {row, col, length, -1, -1}));
    return va_array_length(graph->nodes) - 1;
} // add_node


static inline size_t
add_edge(VA_Allocator const allocator, Graph* const graph, size_t const sink, size_t const next, VA_Variant const variant)
{
    va_array_append(allocator, graph->edges, ((Edge) {sink, next, variant}));
    return va_array_length(graph->edges) - 1;
} // add_edge


static size_t
to_dot(Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    (void) observed;
    size_t count = 0;
    //printf("digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\nsi[shape=point,width=.1]\n");
    //printf("si->s%u\n", graph.source);
    for (size_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        if (graph.nodes[i].edges == (uint32_t) -1)
        {
            // FIXME: always node 0
            //printf("s%zu[label=\"(%u, %u, %u)\",peripheries=2]\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length);
        } // if
        else
        {
            //printf("s%zu[label=\"(%u, %u, %u)\"]\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length);
        } // else

        for (uint32_t j = i; j != (uint32_t) -1; j = graph.nodes[j].lambda)
        {
            for (uint32_t k = graph.nodes[j].edges; k != (uint32_t) -1; k = graph.edges[k].next)
            {
                //printf("s%zu->s%u[label=\"%u:%u/%.*s\"]\n", i, graph.edges[k].sink, graph.edges[k].variant.start, graph.edges[k].variant.end, (int) graph.edges[k].variant.obs_end - graph.edges[k].variant.obs_start, observed + graph.edges[k].variant.obs_start);
                count += 1;
            } // for
        } // for
    } // for
    //printf("}\n");
    return count;
} // to_dot


static void
destroy(VA_Allocator const allocator, Graph* const graph)
{
    graph->nodes = va_array_destroy(allocator, graph->nodes);
    graph->edges = va_array_destroy(allocator, graph->edges);
} // destroy


static size_t
doit(VA_Allocator const allocator, Graph* const graph, VA_LCS_Node const sink, VA_LCS_Node* const sources, size_t const lcs_idx, size_t const shift)
{
    size_t idx = 0;
    for (size_t i = 0; i < va_array_length(sources); ++i)
    {
        if (sources[i].row + sources[i].length < sink.row + sink.length && sources[i].col + sources[i].length < sink.col + sink.length)
        {
            VA_Variant const variant = {sources[i].row + sources[i].length, sink.row + sink.length - 1, sources[i].col + sources[i].length - shift, sink.col + sink.length - 1 - shift};
            //printf("(%u, %u, %u) -> (%u, %u, %u)\n", sources[i].row, sources[i].col, sources[i].length, sink.row, sink.col, sink.length);

            if (sources[i].incoming == lcs_idx)
            {
                //printf("SPLIT: (%u, %u, %u)\n", sources[i].row, sources[i].col, sources[i].length);
                size_t const split_idx = sources[i].idx;
                sources[i].idx = add_node(allocator, graph, sources[i].row, sources[i].col, graph->nodes[split_idx].length);
                sources[i].incoming = 0;

                // lambda-edge
                graph->nodes[sources[i].idx].lambda = split_idx;

                graph->nodes[sources[i].idx].edges = add_edge(allocator, graph, sink.idx, graph->nodes[sources[i].idx].edges , variant);

                graph->nodes[split_idx].row += sources[i].length;
                graph->nodes[split_idx].col += sources[i].length;
                graph->nodes[split_idx].length -= sources[i].length;
            } // if
            else
            {
                if (sources[i].idx == (uint32_t) -1)
                {
                    sources[i].idx = add_node(allocator, graph, sources[i].row, sources[i].col, sources[i].length);
                } // if
                graph->nodes[sources[i].idx].edges = add_edge(allocator, graph, sink.idx, graph->nodes[sources[i].idx].edges, variant);
            } // else
            idx = i + 1;
        } // if
    } // for
    return idx;
} // doit


static Graph
build_graph(VA_Allocator const allocator,
            size_t const len_ref, char const reference[static len_ref],
            size_t const len_obs, char const observed[static len_obs],
            size_t const len_lcs, VA_LCS_Node* lcs_nodes[static len_lcs],
            size_t const shift)
{
    (void) reference;
    (void) observed;

    Graph graph = {
        va_array_init(allocator, 256 * len_lcs, sizeof(*graph.nodes)),
        NULL,
        0,
    };

    if (len_lcs == 0 || lcs_nodes == NULL)
    {
        return graph;  // FIXME: create default graph
    } // if

    //printf("%zu\n", len_lcs);

    VA_LCS_Node sink = lcs_nodes[len_lcs - 1][va_array_length(lcs_nodes[len_lcs - 1]) - 1];
    if (sink.row + sink.length == len_ref + shift && sink.col + sink.length == len_obs + shift)
    {
        va_array_header(lcs_nodes[len_lcs - 1])->length -= 1;  // del
        sink.length += 1;
    } // if
    else
    {
        sink = (VA_LCS_Node) {.row = len_ref + shift, .col = len_obs + shift, .length = 1};
    } // else

    //printf("sink: (%u, %u, %u)\n", sink.row, sink.col, sink.length);

    sink.idx = add_node(allocator, &graph, sink.row, sink.col, sink.length - 1);
    size_t const here = doit(allocator, &graph, sink, lcs_nodes[len_lcs - 1], len_lcs - 1, shift);
    if (sink.length > 1)
    {
        sink.length -= 1;
        sink.incoming = len_lcs - 1;
        va_array_insert(allocator, lcs_nodes[len_lcs - 1], sink, here);
    } // if

    for (size_t i = len_lcs - 1; i >= 1; --i)
    {
        //printf("idx %zu\n", i);
        for (size_t j = 0; j < va_array_length(lcs_nodes[i]); ++j)
        {
            if (lcs_nodes[i][j].idx == (uint32_t) -1)
            {
                continue;
            } // if
            size_t const here = doit(allocator, &graph, lcs_nodes[i][j], lcs_nodes[i - 1], i, shift);
            if (lcs_nodes[i][j].length > 1)
            {
                lcs_nodes[i][j].length -= 1;
                if (here > 0)
                {
                    lcs_nodes[i][j].incoming = i;
                } // if
                va_array_insert(allocator, lcs_nodes[i - 1], lcs_nodes[i][j], here);
            } // if
        } // for
        lcs_nodes[i] = va_array_destroy(allocator, lcs_nodes[i]);
    } // for

    VA_LCS_Node source = lcs_nodes[0][0];

    //printf("source: (%u, %u, %u)\n", source.row, source.col, source.length);
    size_t start = 0;
    if (source.row == shift && source.col == shift)
    {
        // assert source.idx != 0
        start = 1;
    } // if
    else
    {
        source = (VA_LCS_Node) {.row = shift, .col = shift, .length = 0};
        source.idx = add_node(allocator, &graph, shift, shift, 0);
    } // else

    for (size_t i = start; i < va_array_length(lcs_nodes[0]); ++i)
    {
        if (lcs_nodes[0][i].idx == (uint32_t) -1)
        {
            continue;
        } // if

        if (source.row < lcs_nodes[0][i].row + lcs_nodes[0][i].length && source.col < lcs_nodes[0][i].col + lcs_nodes[0][i].length)
        {
            VA_Variant const variant = {source.row, lcs_nodes[0][i].row + lcs_nodes[0][i].length - 1, source.col - shift, lcs_nodes[0][i].col + lcs_nodes[0][i].length - 1 - shift};
            //printf("(%u, %u, %u) -> (%u, %u, %u)\n", source.row, source.col, source.length, lcs_nodes[0][i].row, lcs_nodes[0][i].col, lcs_nodes[0][i].length);
            graph.nodes[source.idx].edges = add_edge(allocator, &graph, lcs_nodes[0][i].idx, graph.nodes[source.idx].edges, variant);
        } // if
    } // for
    lcs_nodes[0] = va_array_destroy(allocator, lcs_nodes[0]);

    lcs_nodes = allocator.alloc(allocator.context, lcs_nodes, len_lcs * sizeof(*lcs_nodes), 0);

    graph.source = source.idx;

    return graph;
} // build_graph


static Graph
reorder(VA_Allocator const allocator, Graph const graph)
{
    Graph new_graph =
    {
        va_array_init(allocator, va_array_length(graph.nodes), sizeof(*new_graph.nodes)),
        va_array_init(allocator, va_array_length(graph.edges), sizeof(*new_graph.edges)),
        graph.source,
    };

    for (size_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        add_node(allocator, &new_graph, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length);

        new_graph.nodes[i].lambda = graph.nodes[i].lambda;
        new_graph.nodes[i].edges = -1;
        for (uint32_t j = graph.nodes[i].edges; j != (uint32_t) -1; j = graph.edges[j].next)
        {
            new_graph.nodes[i].edges = add_edge(allocator, &new_graph, graph.edges[j].sink, new_graph.nodes[i].edges, graph.edges[j].variant);
        } // for
    } // for
    return new_graph;
} // reorder


static size_t
bfs_traversal(Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    (void) observed;
    size_t count = 0;
    VA_Bitset* scheduled = va_bitset_init(va_std_allocator, va_array_length(graph.nodes));
    VA_Queue* queue = va_queue_init(va_std_allocator, va_array_length(graph.nodes));

    //printf("digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\nsi[shape=point,width=.1]\n");
    //printf("si->s%u\n", graph.source);

    va_bitset_set(scheduled, graph.source);
    va_queue_enqueue(queue, graph.source);
    while (!va_queue_is_empty(queue))
    {
        uint32_t const source = va_queue_dequeue(queue);
        if (source == (uint32_t) -1)
        {
            printf("UNDERFLOW\n");
        } // if

        /*
        if (graph.nodes[source].edges == (uint32_t) -1)
        {
            printf("s%u[label=\"(%u, %u, %u)\",peripheries=2]\n", source, graph.nodes[source].row, graph.nodes[source].col, graph.nodes[source].length);
        } // if
        else
        {
            printf("s%u[label=\"(%u, %u, %u)\"]\n", source, graph.nodes[source].row, graph.nodes[source].col, graph.nodes[source].length);
        } // else
        */

        for (uint32_t i = source; i != (uint32_t) -1; i = graph.nodes[i].lambda)
        {
            if (!va_bitset_test(scheduled, i))
            {
                va_bitset_set(scheduled, i);
                if (!va_queue_enqueue(queue, i))
                {
                    printf("OVERFLOW (lambda)\n");
                } // if
            } // if
            for (uint32_t j = graph.nodes[i].edges; j != (uint32_t) -1; j = graph.edges[j].next)
            {
                count += 1;
                //printf("s%u->s%u[label=\"%u:%u/%.*s\"]\n", source, graph.edges[j].sink, graph.edges[j].variant.start, graph.edges[j].variant.end, (int) graph.edges[j].variant.obs_end - graph.edges[j].variant.obs_start, observed + graph.edges[j].variant.obs_start);
                if (!va_bitset_test(scheduled, graph.edges[j].sink))
                {
                    va_bitset_set(scheduled, graph.edges[j].sink);
                    if (!va_queue_enqueue(queue, graph.edges[j].sink))
                    {
                        printf("OVERFLOW\n");
                    } // if
                } // if
            } // for
        } // for

    } // while

    //printf("}\n");

    queue = va_queue_destroy(va_std_allocator, queue);
    scheduled = va_bitset_destroy(va_std_allocator, scheduled);
    return count;
} // bfs_traversal


int
main(int argc, char* argv[static argc + 1])
{

    (void) argv;
    char const* const restrict reference = "ATTCTATCTTCTGTCTACATAAGATGTCATACTAGAGGGCATATCTGCAATGTATACATATTATCTTTTCCAGCATGCATTCAGTTGTGTTGGAATAATTTATGTACACCTTTATAAACGCTGAGCCTCACAAGAGCCATGTGCCACGTATTGTTTTCTTACTACTTTTTGGGATACCTGGCACGTAATAGACACTCATTGAAAGTTTCCTAATGAATGAAGTACAAAGATAAAACAAGTTATAGACTGATTCTTTTGAGCTGTCAAGGTTGTAAATAGACTTTTGCTCAATCAATTCAAATGGTGGCAGGTAGTGGGGGTAGAGGGATTGGTATGAAAAACATAAGCTTTCAGAACTCCTGTGTTTATTTTTAGAATGTCAACTGCTTGAGTGTTTTTAACTCTGTGGTATCTGAACTATCTTCTCTAACTGCAGGTTGGGCTCAGATCTGTGATAGAACAGTTTCCTGGGAAGCTTGACTTTGTCCTTGTGGATGGGGGCTGTGTCCTAAGCCATGGCCACAAGCAGTTGATGTGCTTGGCTAGATCTGTTCTCAGTAAGGCGAAGATCTTGCTGCTTGATGAACCCAGTGCTCATTTGGATCCAGTGTGAGTTTCAGATGTTCTGTTACTTAATAGCACAGTGGGAACAGAATCATTATGCCTGCTTCATGGTGACACATATTTCTATTAGGCTGTCATGTCTGCGTGTGGGGGTCTCCCCCAAGATATGAAATAATTGCCCAGTGGAAATGAGCATAAATGCATATTTCCTTGCTAAGAGTCTTGTGTTTTCTTCCGAAGATAGTTTTTAGTTTCATACAAACTCTTCCCCCTTGTCAACACATGATGAAGCTTTTAAATACATGGGCCTAATCTGATCCTTATGATTTGCCTTTGTATCCCATTTATACCATAAGCATGTTTATAGCCCCAAATAAAGAAGTACTGGTGATTCTACATAATGAAAAATGTACTCATTTATTAAAGTTTCTTTGAAATATTTGTCCTGTTTATTTATGGATACTTAGAGTCTACCCCATGGTTGAAAAGCTGATTGTGGCTAACGCTATATCAACATTATGTGAAAAGAACTTAAAGAAATAAGTAATTTAAAGAGATAATAGAACAATAGACATATTATCAAGGTAAATACAGATCATTACTGTTCTGTGATATTATGTGTGGTATTTTCTTTCTTTTCTAGAACATACCAAATAATTAGAAGAACTCTAAAACAAGCATTTGCTGATTGCACAGTAATTCTCTGTGAACACAGGATAGAAGCAATGCTGGAATGCCAACAATTTTTGGTGAGTCTTTATAACTTTACTTAAGATCTCATTGCCCTTGTAATTCTTGATAACAATCTCACATGTGATAGTTCCTGCAAATTGCAACAATGTACAAGTTCTTTTCAAAAATATGTATCATACAGCCATCCAGCTTTACTCAAAATAGCTGCACAAGTTTTTCACTTTGATCTGAGCCATGTGGTGAGGTTGAAATATAGTAAATCTAAAATGGCAGCATATTACTAAGTTATGTTTATAAATAGGATATATATACTTTTTGAGCCCTTTATTTGGGGACCAAGTCATACAAAATACTCTACTGTTTAAGATTTTAAAAAAGGTCCCTGTGATTCTTTCAATAACTAAATGTCCCATGGATGTGGTCTGGGACAGGCCTAGTTGTCTTACAGTCTGATTTATGGTATTAATGACAAAGTTGAGAGGCACATTTCATTTTT";
    char const* const restrict observed = "ATTCTATCTTCTGTCTACATAAGATGTCATACTAGAGGGCATATCTGCAATGTATACATATTATCTTTTCCAGCATGCATTCAGTTGTGTTGGAATAATTTATGTACACCTTTATAAACGCTGAGCCTCACAAGAGCCATGTGCCACGTATTGTTTTCTTACTACTTTTTGGGATACCTGGCACGTAATAGACACTCATTGAAAGTTTCCTAATGAATGAAGTACAAAGATAAAACAAGTTATAGACTGATTCTTTTGAGCTGTCAAGGTTGTAAATAGACTTTTGCTCAATCAATTCAAATGGTGGCAGGTAGTGGGGGTAGAGGGATTGGTATGAAAAACATAAGCTTTCAGAACTCCTGTGTTTATTTTTAGAATGTCAACTGCTTGAGTGTTTTTAACTCTGTGGTATCTGAACTATCTTCTCTAACTGCAGGTGAGTCTTTATAACTTTACTTAAGATCTCATTGCCCTTGTAATTCTTGATAACAATCTCACATGTGATAGTTCCTGCAAATTGCAACAATGTACAAGTTCTTTTCAAAAATATGTATCATACAGCCATCCAGCTTTACTCAAAATAGCTGCACAAGTTTTTCACTTTGATCTGAGCCATGTGGTGAGGTTGAAATATAGTAAATCTAAAATGGCAGCATATTACTAAGTTATGTTTATAAATAGGATATATATACTTTTTGAGCCCTTTATTTGGGGACCAAGTCATACAAAATACTCTACTGTTTAAGATTTTAAAAAAGGTCCCTGTGATTCTTTCAATAACTAAATGTCCCATGGATGTGGTCTGGGACAGGCCTAGTTGTCTTACAGTCTGATTTATGGTATTAATGACAAAGTTGAGAGGCACATTTCATTTTT";

/*
    if (argc < 3)
    {
        fprintf(stderr, "usage: %s reference observed\n", argv[0]);
        return EXIT_FAILURE;
    } // if
    char const* const restrict reference = argv[1];
    char const* const restrict observed = argv[2];
*/
    size_t const len_ref = strlen(reference);
    size_t const len_obs = strlen(observed);

    VA_LCS_Node** lcs_nodes = NULL;
    size_t const len_lcs = va_edit(va_std_allocator, len_ref, reference, len_obs, observed, &lcs_nodes);
    Graph graph = build_graph(va_std_allocator, len_ref, reference, len_obs, observed, len_lcs, lcs_nodes, 0);

    printf("%zu\n", bfs_traversal(graph, len_obs, observed));
    printf("%zu\n", to_dot(graph, len_obs, observed));

    destroy(va_std_allocator, &graph);

    return EXIT_SUCCESS;
} // main
