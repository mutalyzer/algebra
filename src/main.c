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
    uint32_t tail;
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
add_edge(VA_Allocator const allocator, Graph* const graph, size_t const tail, size_t const next, VA_Variant const variant)
{
    va_array_append(allocator, graph->edges, ((Edge) {tail, next, variant}));
    return va_array_length(graph->edges) - 1;
} // add_edge


static size_t
to_dot(Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    (void) observed;
    size_t count = 0;
    printf("digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\nsi[shape=point,width=.1]\n");
    printf("si->s%u\n", graph.source);
    for (size_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        if (graph.nodes[i].edges == (uint32_t) -1)
        {
            // FIXME: always node 0
            printf("s%zu[label=\"(%u, %u, %u)\",peripheries=2]\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length);
        } // if
        else
        {
            printf("s%zu[label=\"(%u, %u, %u)\"]\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length);
        } // else

        for (uint32_t j = i; j != (uint32_t) -1; j = graph.nodes[j].lambda)
        {
            for (uint32_t k = graph.nodes[j].edges; k != (uint32_t) -1; k = graph.edges[k].next)
            {
                printf("s%zu->s%u[label=\"%u:%u/%.*s\"]\n", i, graph.edges[k].tail, graph.edges[k].variant.start, graph.edges[k].variant.end, (int) graph.edges[k].variant.obs_end - graph.edges[k].variant.obs_start, observed + graph.edges[k].variant.obs_start);
                count += 1;
            } // for
        } // for
    } // for
    printf("}\n");
    return count;
} // to_dot


static void
destroy(VA_Allocator const allocator, Graph* const graph)
{
    graph->nodes = va_array_destroy(allocator, graph->nodes);
    graph->edges = va_array_destroy(allocator, graph->edges);
} // destroy


static size_t
doit(VA_Allocator const allocator, Graph* const graph, VA_LCS_Node const tail, VA_LCS_Node* const heads, size_t const lcs_idx, size_t const shift)
{
    size_t idx = 0;
    for (size_t i = 0; i < va_array_length(heads); ++i)
    {
        if (heads[i].row + heads[i].length < tail.row + tail.length && heads[i].col + heads[i].length < tail.col + tail.length)
        {
            VA_Variant const variant = {heads[i].row + heads[i].length, tail.row + tail.length - 1, heads[i].col + heads[i].length - shift, tail.col + tail.length - 1 - shift};
            //printf("(%u, %u, %u) -> (%u, %u, %u)\n", heads[i].row, heads[i].col, heads[i].length, tail.row, tail.col, tail.length);

            if (heads[i].incoming == lcs_idx)
            {
                //printf("SPLIT: (%u, %u, %u)\n", heads[i].row, heads[i].col, heads[i].length);
                size_t const split_idx = heads[i].idx;
                heads[i].idx = add_node(allocator, graph, heads[i].row, heads[i].col, graph->nodes[split_idx].length);
                heads[i].incoming = 0;

                // lambda-edge
                graph->nodes[heads[i].idx].lambda = split_idx;

                graph->nodes[heads[i].idx].edges = add_edge(allocator, graph, tail.idx, graph->nodes[heads[i].idx].edges , variant);

                graph->nodes[split_idx].row += heads[i].length;
                graph->nodes[split_idx].col += heads[i].length;
                graph->nodes[split_idx].length -= heads[i].length;
            } // if
            else
            {
                if (heads[i].idx == (uint32_t) -1)
                {
                    heads[i].idx = add_node(allocator, graph, heads[i].row, heads[i].col, heads[i].length);
                } // if
                graph->nodes[heads[i].idx].edges = add_edge(allocator, graph, tail.idx, graph->nodes[heads[i].idx].edges, variant);
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
        //va_array_init(allocator, 256 * len_lcs, sizeof(*graph.nodes)),  // FIXME: we know the number of lcs nodes
        NULL,
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
            new_graph.nodes[i].edges = add_edge(allocator, &new_graph, graph.edges[j].tail, new_graph.nodes[i].edges, graph.edges[j].variant);
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

    /*
    printf("digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\nsi[shape=point,width=.1]\n");
    printf("si->s%u\n", graph.source);
    */

    va_bitset_set(scheduled, graph.source);
    va_queue_enqueue(queue, graph.source);
    while (!va_queue_is_empty(queue))
    {
        uint32_t const head = va_queue_dequeue(queue);

        /*
        if (graph.nodes[head].edges == (uint32_t) -1)
        {
            printf("s%u[label=\"(%u, %u, %u)\",peripheries=2]\n", head, graph.nodes[head].row, graph.nodes[head].col, graph.nodes[head].length);
        } // if
        else
        {
            printf("s%u[label=\"(%u, %u, %u)\"]\n", head, graph.nodes[head].row, graph.nodes[head].col, graph.nodes[head].length);
        } // else
        */

        for (uint32_t i = head; i != (uint32_t) -1; i = graph.nodes[i].lambda)
        {
            if (!va_bitset_test(scheduled, i))
            {
                va_bitset_set(scheduled, i);
                va_queue_enqueue(queue, i);  // OVERFLOW
            } // if
            for (uint32_t j = graph.nodes[i].edges; j != (uint32_t) -1; j = graph.edges[j].next)
            {
                count += 1;
                //printf("s%u->s%u[label=\"%u:%u/%.*s\"]\n", head, graph.edges[j].tail, graph.edges[j].variant.start, graph.edges[j].variant.end, (int) graph.edges[j].variant.obs_end - graph.edges[j].variant.obs_start, observed + graph.edges[j].variant.obs_start);
                if (!va_bitset_test(scheduled, graph.edges[j].tail))
                {
                    va_bitset_set(scheduled, graph.edges[j].tail);
                    va_queue_enqueue(queue, graph.edges[j].tail);  // OVERFLOW
                } // if
            } // for
        } // for
    } // while

    //printf("}\n");

    queue = va_queue_destroy(va_std_allocator, queue);
    scheduled = va_bitset_destroy(va_std_allocator, scheduled);
    return count;
} // bfs_traversal


static size_t
bfs_traversal_alt(VA_Allocator const allocator, Graph const graph)
{
    struct
    {
        uint32_t depth;
        uint32_t next;
    }* queue = allocator.alloc(allocator.context, NULL, 0, sizeof(*queue) * va_array_length(graph.nodes));

    for (size_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        queue[i].depth = 0;
        queue[i].next = -1;
    } // for
    uint32_t head = graph.source;
    uint32_t tail = graph.source;

    size_t count = 0;

    while (head != (uint32_t) -1)
    {
        printf("pop %u\n", head);
        for (uint32_t i = head; i != (uint32_t) -1; i = graph.nodes[i].lambda)
        {
            if (i != head)
            {
                printf("lambda %u\n", i);
            } // if
            for (uint32_t j = graph.nodes[i].edges; j != (uint32_t) -1; j = graph.edges[j].next)
            {
                count += 1;

                if (queue[graph.edges[j].tail].depth > 0)
                {
                    printf("skip %u\n", graph.edges[j].tail);
                    continue;
                } // if

                printf("push %u\n", graph.edges[j].tail);
                queue[graph.edges[j].tail].depth = queue[i].depth + 1;
                queue[tail].next = graph.edges[j].tail;
                tail = graph.edges[j].tail;
            } // for
        } // for
        head = queue[head].next;

        printf("  # \tdepth\tnext\n");
        for (size_t i = 0; i < va_array_length(graph.nodes); ++i)
        {
            printf("%3zu:\t%5u\t%4d\n", i, queue[i].depth, queue[i].next);
        } // for
        printf("head: %d\ntail: %d\n", head, tail);
    } // while

    queue = allocator.alloc(allocator.context, queue, sizeof(*queue) * va_array_length(graph.nodes), 0);

    return count;
} // bfs_traversal_alt


int
main(int argc, char* argv[static argc + 1])
{
/*
    (void) argv;
    char const* const restrict reference = "ATTCTATCTTCTGTCTACATAAGATGTCATACTAGAGGGCATATCTGCAATGTATACATATTATCTTTTCCAGCATGCATTCAGTTGTGTTGGAATAATTTATGTACACCTTTATAAACGCTGAGCCTCACAAGAGCCATGTGCCACGTATTGTTTTCTTACTACTTTTTGGGATACCTGGCACGTAATAGACACTCATTGAAAGTTTCCTAATGAATGAAGTACAAAGATAAAACAAGTTATAGACTGATTCTTTTGAGCTGTCAAGGTTGTAAATAGACTTTTGCTCAATCAATTCAAATGGTGGCAGGTAGTGGGGGTAGAGGGATTGGTATGAAAAACATAAGCTTTCAGAACTCCTGTGTTTATTTTTAGAATGTCAACTGCTTGAGTGTTTTTAACTCTGTGGTATCTGAACTATCTTCTCTAACTGCAGGTTGGGCTCAGATCTGTGATAGAACAGTTTCCTGGGAAGCTTGACTTTGTCCTTGTGGATGGGGGCTGTGTCCTAAGCCATGGCCACAAGCAGTTGATGTGCTTGGCTAGATCTGTTCTCAGTAAGGCGAAGATCTTGCTGCTTGATGAACCCAGTGCTCATTTGGATCCAGTGTGAGTTTCAGATGTTCTGTTACTTAATAGCACAGTGGGAACAGAATCATTATGCCTGCTTCATGGTGACACATATTTCTATTAGGCTGTCATGTCTGCGTGTGGGGGTCTCCCCCAAGATATGAAATAATTGCCCAGTGGAAATGAGCATAAATGCATATTTCCTTGCTAAGAGTCTTGTGTTTTCTTCCGAAGATAGTTTTTAGTTTCATACAAACTCTTCCCCCTTGTCAACACATGATGAAGCTTTTAAATACATGGGCCTAATCTGATCCTTATGATTTGCCTTTGTATCCCATTTATACCATAAGCATGTTTATAGCCCCAAATAAAGAAGTACTGGTGATTCTACATAATGAAAAATGTACTCATTTATTAAAGTTTCTTTGAAATATTTGTCCTGTTTATTTATGGATACTTAGAGTCTACCCCATGGTTGAAAAGCTGATTGTGGCTAACGCTATATCAACATTATGTGAAAAGAACTTAAAGAAATAAGTAATTTAAAGAGATAATAGAACAATAGACATATTATCAAGGTAAATACAGATCATTACTGTTCTGTGATATTATGTGTGGTATTTTCTTTCTTTTCTAGAACATACCAAATAATTAGAAGAACTCTAAAACAAGCATTTGCTGATTGCACAGTAATTCTCTGTGAACACAGGATAGAAGCAATGCTGGAATGCCAACAATTTTTGGTGAGTCTTTATAACTTTACTTAAGATCTCATTGCCCTTGTAATTCTTGATAACAATCTCACATGTGATAGTTCCTGCAAATTGCAACAATGTACAAGTTCTTTTCAAAAATATGTATCATACAGCCATCCAGCTTTACTCAAAATAGCTGCACAAGTTTTTCACTTTGATCTGAGCCATGTGGTGAGGTTGAAATATAGTAAATCTAAAATGGCAGCATATTACTAAGTTATGTTTATAAATAGGATATATATACTTTTTGAGCCCTTTATTTGGGGACCAAGTCATACAAAATACTCTACTGTTTAAGATTTTAAAAAAGGTCCCTGTGATTCTTTCAATAACTAAATGTCCCATGGATGTGGTCTGGGACAGGCCTAGTTGTCTTACAGTCTGATTTATGGTATTAATGACAAAGTTGAGAGGCACATTTCATTTTT";
    char const* const restrict observed = "ATTCTATCTTCTGTCTACATAAGATGTCATACTAGAGGGCATATCTGCAATGTATACATATTATCTTTTCCAGCATGCATTCAGTTGTGTTGGAATAATTTATGTACACCTTTATAAACGCTGAGCCTCACAAGAGCCATGTGCCACGTATTGTTTTCTTACTACTTTTTGGGATACCTGGCACGTAATAGACACTCATTGAAAGTTTCCTAATGAATGAAGTACAAAGATAAAACAAGTTATAGACTGATTCTTTTGAGCTGTCAAGGTTGTAAATAGACTTTTGCTCAATCAATTCAAATGGTGGCAGGTAGTGGGGGTAGAGGGATTGGTATGAAAAACATAAGCTTTCAGAACTCCTGTGTTTATTTTTAGAATGTCAACTGCTTGAGTGTTTTTAACTCTGTGGTATCTGAACTATCTTCTCTAACTGCAGGTGAGTCTTTATAACTTTACTTAAGATCTCATTGCCCTTGTAATTCTTGATAACAATCTCACATGTGATAGTTCCTGCAAATTGCAACAATGTACAAGTTCTTTTCAAAAATATGTATCATACAGCCATCCAGCTTTACTCAAAATAGCTGCACAAGTTTTTCACTTTGATCTGAGCCATGTGGTGAGGTTGAAATATAGTAAATCTAAAATGGCAGCATATTACTAAGTTATGTTTATAAATAGGATATATATACTTTTTGAGCCCTTTATTTGGGGACCAAGTCATACAAAATACTCTACTGTTTAAGATTTTAAAAAAGGTCCCTGTGATTCTTTCAATAACTAAATGTCCCATGGATGTGGTCTGGGACAGGCCTAGTTGTCTTACAGTCTGATTTATGGTATTAATGACAAAGTTGAGAGGCACATTTCATTTTT";
*/

    if (argc < 3)
    {
        fprintf(stderr, "usage: %s reference observed\n", argv[0]);
        return EXIT_FAILURE;
    } // if
    char const* const restrict reference = argv[1];
    char const* const restrict observed = argv[2];

    size_t const len_ref = strlen(reference);
    size_t const len_obs = strlen(observed);

    VA_LCS_Node** lcs_nodes = NULL;
    size_t const len_lcs = va_edit(va_std_allocator, len_ref, reference, len_obs, observed, &lcs_nodes);
    /*
    printf("%zu\n", len_lcs);
    for (size_t i = 0; i < len_lcs; ++i)
    {
        printf("%zu:  ", i);
        for (size_t j = 0; j < va_array_length(lcs_nodes[i]); ++j)
        {
            printf("(%u, %u, %u) ", lcs_nodes[i][j].row, lcs_nodes[i][j].col, lcs_nodes[i][j].length);
        } // for
        printf("\n");
    } // for
    */

    Graph graph = build_graph(va_std_allocator, len_ref, reference, len_obs, observed, len_lcs, lcs_nodes, 0);

    printf("%zu\n", to_dot(graph, len_obs, observed));

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

    printf("%zu\n", bfs_traversal_alt(va_std_allocator, graph));

    destroy(va_std_allocator, &graph);

    return EXIT_SUCCESS;
} // main
