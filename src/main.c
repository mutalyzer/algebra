#include <stdio.h>      // stderr, fprintf, printf
#include <stdlib.h>     // EXIT_*
#include <string.h>     // strlen

#include "../include/alloc.h"           // VA_Allocator
#include "../include/array.h"           // va_array_*
#include "../include/edit.h"            // va_edit
#include "../include/std_alloc.h"       // va_std_allocator
#include "../include/variant.h"         // VA_Variant


typedef struct
{
    size_t node;
    VA_Variant variant;
} Edge;


typedef struct
{
    size_t row;
    size_t col;
    size_t length;
    Edge* edges;
} Node;


typedef struct
{
    Node* nodes;
} Graph;


static inline size_t
add_node(VA_Allocator const allocator, Graph* const graph, size_t const row, size_t const col, size_t const length)
{
    va_array_append(allocator, graph->nodes, ((Node) {row, col, length, NULL}));
    return va_array_length(graph->nodes);
} // add_node


static inline void
add_edge(VA_Allocator const allocator, Graph* const graph, size_t const source, size_t const sink, VA_Variant const variant)
{
    va_array_append(allocator, graph->nodes[source].edges, ((Edge) {sink, variant}));
} // add_edge


static void
to_dot(Graph const* const graph, char const observed[static 1])
{
    if (graph->nodes == NULL)
    {
        return;
    } // if
    for (size_t i = 0; i < va_array_length(graph->nodes); ++i)
    {
        printf("s%zu[label=\"(%zu, %zu, %zu)\"]\n", i, graph->nodes[i].row, graph->nodes[i].col, graph->nodes[i].length);
    } // for

    for (size_t i = 0; i < va_array_length(graph->nodes); ++i)
    {
        for (size_t j = 0; j < va_array_length(graph->nodes[i].edges); ++j)
        {
            size_t const len = graph->nodes[i].edges[j].variant.obs_end - graph->nodes[i].edges[j].variant.obs_start;
            printf("s%zu->s%zu[label=\"%u:%u/%.*s\"]\n", i, graph->nodes[i].edges[j].node, graph->nodes[i].edges[j].variant.start, graph->nodes[i].edges[j].variant.end, (int) len, observed + graph->nodes[i].edges[j].variant.obs_start);
        } // for
    } // for
} // to_dot


static void
destroy(VA_Allocator const allocator, Graph* const graph)
{
    if (graph->nodes == NULL)
    {
        return;
    } // if
    for (size_t i = 0; i < va_array_length(graph->nodes); ++i)
    {
        graph->nodes[i].edges = va_array_destroy(allocator, graph->nodes[i].edges);
    } // for
    graph->nodes = va_array_destroy(allocator, graph->nodes);
} // destroy


static void
doit(VA_Allocator const allocator, Graph* const graph, VA_LCS_Node const sink, VA_LCS_Node* const sources, size_t const shift)
{
    for (size_t i = 0; i < va_array_length(sources); ++i)
    {
        if (sources[i].row + sources[i].length < sink.row + sink.length && sources[i].col + sources[i].length < sink.col + sink.length)
        {
            VA_Variant const variant = {sources[i].row + sources[i].length, sink.row + sink.length - 1, sources[i].col + sources[i].length - shift, sink.col + sink.length - 1 - shift};
            printf("(%zu, %zu, %zu) -> (%zu, %zu, %zu)\n", sources[i].row, sources[i].col, sources[i].length, sink.row, sink.col, sink.length);

            if (sources[i].idx == 0)
            {
                sources[i].idx = add_node(allocator, graph, sources[i].row, sources[i].col, sources[i].length);
            } // if

            add_edge(allocator, graph, sources[i].idx - 1, sink.idx - 1, variant);

        } // if
    } // for
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

    Graph graph = {NULL};

    if (len_lcs == 0 || lcs_nodes == NULL)
    {
        return graph;  // FIXME: create default graph
    } // if

    printf("%zu\n", len_lcs);

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

    printf("sink: (%zu, %zu, %zu)\n", sink.row, sink.col, sink.length);

    sink.idx = add_node(allocator, &graph, sink.row, sink.col, sink.length - 1);
    doit(allocator, &graph, sink, lcs_nodes[len_lcs - 1], shift);

    for (size_t i = 0; i < len_lcs; ++i)
    {
        size_t const idx = len_lcs - 1 - i;
        printf("%zu (%zu): ", idx, va_array_length(lcs_nodes[idx]));
        for (size_t j = 0; j < va_array_length(lcs_nodes[idx]); ++j)
        {
            printf("(%zu, %zu, %zu) ", lcs_nodes[idx][j].row, lcs_nodes[idx][j].col, lcs_nodes[idx][j].length);
        } // for
        printf("\n");
        lcs_nodes[idx] = va_array_destroy(allocator, lcs_nodes[idx]);
    } // for

    lcs_nodes = allocator.alloc(allocator.context, lcs_nodes, len_lcs * sizeof(*lcs_nodes), 0);

    return graph;
} // build_graph


int
main(int argc, char* argv[static argc + 1])
{
    if (argc < 3)
    {
        fprintf(stderr, "usage: %s reference observed\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    size_t const m = strlen(argv[1]);
    size_t const n = strlen(argv[2]);

    VA_LCS_Node** lcs_nodes = NULL;
    size_t const len_lcs = va_edit(va_std_allocator, m, argv[1], n, argv[2], &lcs_nodes);
    Graph graph = build_graph(va_std_allocator, m, argv[1], n, argv[2], len_lcs, lcs_nodes, 0);

    to_dot(&graph, argv[2]);

    destroy(va_std_allocator, &graph);

    return EXIT_SUCCESS;
} // main
