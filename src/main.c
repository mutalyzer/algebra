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
    printf("digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\nsi[shape=point,width=.1]\n");
    if (graph->nodes == NULL)
    {
        printf("}\n");
        return;
    } // if
    for (size_t i = 0; i < va_array_length(graph->nodes); ++i)
    {
        if (va_array_length(graph->nodes[i].edges) == 0)
        {
            printf("s%zu[label=\"(%zu, %zu, %zu)\",peripheries=2]\n", i, graph->nodes[i].row, graph->nodes[i].col, graph->nodes[i].length);
        } // if
        else
        {
            printf("s%zu[label=\"(%zu, %zu, %zu)\"]\n", i, graph->nodes[i].row, graph->nodes[i].col, graph->nodes[i].length);
        } // else
    } // for

    for (size_t i = 0; i < va_array_length(graph->nodes); ++i)
    {
        for (size_t j = 0; j < va_array_length(graph->nodes[i].edges); ++j)
        {
            size_t const len = graph->nodes[i].edges[j].variant.obs_end - graph->nodes[i].edges[j].variant.obs_start;
            printf("s%zu->s%zu[label=\"%u:%u/%.*s\"]\n", i, graph->nodes[i].edges[j].node, graph->nodes[i].edges[j].variant.start, graph->nodes[i].edges[j].variant.end, (int) len, observed + graph->nodes[i].edges[j].variant.obs_start);
        } // for
    } // for
    printf("}\n");
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


static size_t
doit(VA_Allocator const allocator, Graph* const graph, VA_LCS_Node const sink, VA_LCS_Node* const sources, size_t const lcs_idx, size_t const shift)
{
    size_t idx = 0;
    for (size_t i = 0; i < va_array_length(sources); ++i)
    {
        if (sources[i].row + sources[i].length < sink.row + sink.length && sources[i].col + sources[i].length < sink.col + sink.length)
        {
            VA_Variant const variant = {sources[i].row + sources[i].length, sink.row + sink.length - 1, sources[i].col + sources[i].length - shift, sink.col + sink.length - 1 - shift};
            printf("(%zu, %zu, %zu) -> (%zu, %zu, %zu)\n", sources[i].row, sources[i].col, sources[i].length, sink.row, sink.col, sink.length);

            if (sources[i].incoming == lcs_idx)
            {
                printf("SPLIT: (%zu, %zu, %zu)\n", sources[i].row, sources[i].col, sources[i].length);
                size_t const split_idx = sources[i].idx - 1;
                sources[i].idx = add_node(allocator, graph, sources[i].row, sources[i].col, graph->nodes[split_idx].length);

                graph->nodes[sources[i].idx - 1].edges = va_array_init(allocator, va_array_length(graph->nodes[split_idx].edges) + 1, sizeof(*graph->nodes[sources[i].idx].edges));
                for (size_t j = 0; j < va_array_length(graph->nodes[split_idx].edges); ++j)
                {
                    va_array_append(allocator, graph->nodes[sources[i].idx - 1].edges, graph->nodes[split_idx].edges[j]);
                } // for
                add_edge(allocator, graph, sources[i].idx - 1, sink.idx - 1, variant);

                graph->nodes[split_idx].row += sources[i].length;
                graph->nodes[split_idx].col += sources[i].length;
                graph->nodes[split_idx].length -= sources[i].length;
            } // if
            else
            {
                if (sources[i].idx == 0)
                {
                    sources[i].idx = add_node(allocator, graph, sources[i].row, sources[i].col, sources[i].length);
                } // if
                add_edge(allocator, graph, sources[i].idx - 1, sink.idx - 1, variant);
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
    size_t const here = doit(allocator, &graph, sink, lcs_nodes[len_lcs - 1], len_lcs - 1, shift);
    if (sink.length > 1)
    {
        sink.length -= 1;
        sink.incoming = len_lcs - 1;
        va_array_insert(allocator, lcs_nodes[len_lcs - 1], sink, here);
    } // if

    for (size_t i = len_lcs - 1; i >= 1; --i)
    {
        printf("idx %zu\n", i);
        for (size_t j = 0; j < va_array_length(lcs_nodes[i]); ++j)
        {
            if (lcs_nodes[i][j].idx == 0)
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

    printf("source: (%zu, %zu, %zu)\n", source.row, source.col, source.length);
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
        if (lcs_nodes[0][i].idx == 0)
        {
            continue;
        } // if

        if (source.row < lcs_nodes[0][i].row + lcs_nodes[0][i].length && source.col < lcs_nodes[0][i].col + lcs_nodes[0][i].length)
        {
            VA_Variant const variant = {source.row, lcs_nodes[0][i].row + lcs_nodes[0][i].length - 1, source.col - shift, lcs_nodes[0][i].col + lcs_nodes[0][i].length - 1 - shift};
            printf("(%zu, %zu, %zu) -> (%zu, %zu, %zu)\n", source.row, source.col, source.length, lcs_nodes[0][i].row, lcs_nodes[0][i].col, lcs_nodes[0][i].length);
            add_edge(allocator, &graph, source.idx - 1, lcs_nodes[0][i].idx - 1, variant);
        } // if
    } // for
    lcs_nodes[0] = va_array_destroy(allocator, lcs_nodes[0]);

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
