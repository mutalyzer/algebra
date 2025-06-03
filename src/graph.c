#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // intmax_t, uint32_t
#include <stdio.h>      // stderr, fprintf, printf
#include <stdlib.h>     // EXIT_*, atoi, rand, srand
#include <string.h>     // strlen, strncmp

#include "../include/alloc.h"           // VA_Allocator
#include "../include/array.h"           // va_array_*
#include "../include/edit.h"            // VA_LCS_Node, va_edit
#include "../include/extract.h"         // canonical, local_supremal
#include "../include/std_alloc.h"       // va_std_allocator
#include "../include/variant.h"         // VA_Variant
// #include "../include/graph.h"           // Edge, Graph, Node


#define MAX(lhs, rhs) (((lhs) > (rhs)) ? (lhs) : (rhs))
#define MIN(lhs, rhs) (((lhs) < (rhs)) ? (lhs) : (rhs))

#define print_variant(variant, observed) variant.start, variant.end, (int) variant.obs_end - variant.obs_start, observed + variant.obs_start
#define VAR_FMT "%u:%u/%.*s"


static uint32_t const GVA_NULL = UINT32_MAX;


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


inline void
destroy(VA_Allocator const allocator, Graph* const graph)
{
    graph->nodes = va_array_destroy(allocator, graph->nodes);
    graph->edges = va_array_destroy(allocator, graph->edges);
} // destroy


Graph
build_graph(VA_Allocator const allocator,
            size_t const len_ref, size_t const len_obs,
            size_t const len_lcs, VA_LCS_Node* lcs_nodes[static len_lcs],
            size_t const shift, bool const debug)
{
    Graph graph = {.nodes = NULL, .edges = NULL};

    if (len_lcs == 0 || lcs_nodes == NULL)
    {
        graph.source = add_node(allocator, &graph, shift, shift, 0);
        if (len_ref == 0 && len_obs == 0)
        {
            graph.supremal = (VA_Variant) {0, 0, 0, 0};
            return graph;
        } // if
        uint32_t const sink = add_node(allocator, &graph, len_ref, len_obs, 0);
        graph.supremal = (VA_Variant) {shift, shift + len_ref, 0, len_obs};
        graph.nodes[graph.source].edges = add_edge(allocator, &graph, sink, -1, graph.supremal);
        return graph;
    } // if

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

    uint32_t max_sink = 0;
    sink.idx = add_node(allocator, &graph, sink.row, sink.col, sink.length - 1);
    if (debug)
        printf("MAKE NODE %u: (%u, %u, %u)\n", sink.idx, sink.row, sink.col, sink.length - 1);
    size_t here = 0;
    VA_LCS_Node* const heads = lcs_nodes[len_lcs - 1];
    for (size_t i = 0; i < va_array_length(heads); ++i)
    {
        eq_count += 1;
        fprintf(stderr, "(%u, %u, %u) vs (%u, %u, %u)\n", heads[i].row, heads[i].col, heads[i].length, sink.row, sink.col, sink.length);

        if (heads[i].row + heads[i].length < sink.row + sink.length && heads[i].col + heads[i].length < sink.col + sink.length)
        {
            VA_Variant const variant = {heads[i].row + heads[i].length, sink.row + sink.length - 1, heads[i].col + heads[i].length - shift, sink.col + sink.length - 1 - shift};

            max_sink = MAX(max_sink, sink.row + sink.length - 1);
            heads[i].idx = add_node(allocator, &graph, heads[i].row, heads[i].col, heads[i].length);
            if (debug)
                printf("MAKE NODE %u: (%u, %u, %u)\n", heads[i].idx, heads[i].row, heads[i].col, heads[i].length);
            graph.nodes[heads[i].idx].edges = add_edge(allocator, &graph, sink.idx, graph.nodes[heads[i].idx].edges, variant);
            if (debug)
                printf("MAKE EDGE (%u %u %u) -> (%u %u %u): " VAR_FMT "\n", heads[i].row, heads[i].col, heads[i].length, sink.row, sink.col, sink.length, print_variant(variant, "0123456789"));

            here = i + 1;
        } // if
    } // for

    if (sink.length > 1)
    {
        sink.length -= 1;
        sink.incoming = len_lcs - 1;
        va_array_insert(allocator, lcs_nodes[len_lcs - 1], sink, here);
        if (debug)
            printf("INSERT (%u %u %u) here: %zu\n", sink.row, sink.col, sink.length, here);
    } // if

    for (size_t i = len_lcs - 1; i >= 1; --i)
    {
        for (size_t j = 0; j < va_array_length(lcs_nodes[i]); ++j)
        {
            if (lcs_nodes[i][j].idx == GVA_NULL)
            {
                continue;
            } // if

            here = 0;
            VA_LCS_Node const tail = lcs_nodes[i][j];
            VA_LCS_Node* const heads = lcs_nodes[i - 1];
            for (size_t k = 0; k < va_array_length(heads); ++k)
            {
                /*
                if (heads[k].idx != GVA_NULL && graph.nodes[heads[k].idx].length != heads[k].length &&
                    tail.length != graph.nodes[tail.idx].length)
                {
                    continue;
                } // if
                */
                eq_count += 1;
                fprintf(stderr, "(%u, %u, %u) vs (%u, %u, %u)\n", heads[k].row, heads[k].col, heads[k].length, tail.row, tail.col, tail.length);

                if (heads[k].row + heads[k].length < tail.row + tail.length &&
                    heads[k].col + heads[k].length < tail.col + tail.length)
                {
                    max_sink = MAX(max_sink, tail.row + tail.length - 1);

                    VA_Variant const variant = {heads[k].row + heads[k].length, tail.row + tail.length - 1, heads[k].col + heads[k].length - shift, tail.col + tail.length - 1 - shift};

                    if (heads[k].incoming == i)
                    {
                        if (debug)
                            fprintf(stderr, "SPLIT at %zu\n", i);
                        uint32_t const split_idx = heads[k].idx;
                        heads[k].idx = add_node(allocator, &graph, heads[k].row, heads[k].col, heads[k].length);
                        if (debug)
                            printf("MAKE NODE %u: (%u, %u, %u)\n", heads[k].idx, heads[k].row, heads[k].col, heads[k].length);
                        heads[k].incoming = 0;

                        // lambda-edge
                        graph.nodes[heads[k].idx].lambda = split_idx;

                        graph.nodes[split_idx].row += heads[k].length;
                        graph.nodes[split_idx].col += heads[k].length;
                        graph.nodes[split_idx].length -= heads[k].length;
                        if (debug)
                            printf("RELABEL NODE %u: (%u, %u, %u)\n", split_idx, graph.nodes[split_idx].row,graph.nodes[split_idx].col, graph.nodes[split_idx].length);
                    } // if
                    else if (heads[k].idx == GVA_NULL)
                    {
                        heads[k].idx = add_node(allocator, &graph, heads[k].row, heads[k].col, heads[k].length);
                        if (debug)
                            printf("MAKE NODE %u: (%u, %u, %u)\n", heads[k].idx, heads[k].row, heads[k].col, heads[k].length);
                    } // if

                    graph.nodes[heads[k].idx].edges = add_edge(allocator, &graph, tail.idx, graph.nodes[heads[k].idx].edges, variant);
                    if (debug)
                        printf("MAKE EDGE (%u %u %u) -> (%u %u %u): " VAR_FMT "\n", heads[k].row, heads[k].col, heads[k].length, tail.row, tail.col, tail.length, print_variant(variant, "0123456789"));

                    here = k + 1;
                } // if
            } // for

            if (lcs_nodes[i][j].length > 1)
            {
                lcs_nodes[i][j].length -= 1;
                if (here > 0)
                {
                    lcs_nodes[i][j].incoming = i;
                } // if
                va_array_insert(allocator, lcs_nodes[i - 1], lcs_nodes[i][j], here);
                if (debug)
                    printf("INSERT (%u %u %u) here: %zu\n", lcs_nodes[i][j].row, lcs_nodes[i][j].col, lcs_nodes[i][j].length, here);
            } // if
        } // for
        lcs_nodes[i] = va_array_destroy(allocator, lcs_nodes[i]);
    } // for

    VA_LCS_Node source = lcs_nodes[0][0];

    size_t start = 0;
    if (source.row == shift && source.col == shift)
    {
        start = 1;
    } // if
    else
    {
        source = (VA_LCS_Node) {.row = shift, .col = shift, .length = 0};
        source.idx = add_node(allocator, &graph, source.row, source.col, source.length);
        if (debug)
            printf("MAKE NODE %u: (%u, %u, %u)\n", source.idx, source.row, source.col, source.length);
    } // else

    for (size_t i = start; i < va_array_length(lcs_nodes[0]); ++i)
    {
        if (lcs_nodes[0][i].idx == GVA_NULL)
        {
            continue;
        } // if

        eq_count += 1;
        fprintf(stderr, "(%u, %u, %u) vs (%u, %u, %u)\n", source.row, source.col, source.length, lcs_nodes[0][i].row, lcs_nodes[0][i].col, lcs_nodes[0][i].length);

        if (source.row < lcs_nodes[0][i].row + lcs_nodes[0][i].length && source.col < lcs_nodes[0][i].col + lcs_nodes[0][i].length)
        {
            max_sink = MAX(max_sink, lcs_nodes[0][i].row + lcs_nodes[0][i].length - 1);

            VA_Variant const variant = {source.row, lcs_nodes[0][i].row + lcs_nodes[0][i].length - 1, source.col - shift, lcs_nodes[0][i].col + lcs_nodes[0][i].length - 1 - shift};

            graph.nodes[source.idx].edges = add_edge(allocator, &graph, lcs_nodes[0][i].idx, graph.nodes[source.idx].edges, variant);
            if (debug)
                printf("MAKE EDGE (%u %u %u) -> (%u %u %u): " VAR_FMT "\n", source.row, source.col, source.length, lcs_nodes[0][i].row, lcs_nodes[0][i].col, lcs_nodes[0][i].length, print_variant(variant, "0123456789"));
        } // if
    } // for
    lcs_nodes[0] = va_array_destroy(allocator, lcs_nodes[0]);
    lcs_nodes = allocator.alloc(allocator.context, lcs_nodes, len_lcs * sizeof(*lcs_nodes), 0);

    uint32_t min_source = max_sink;
    for (uint32_t i = graph.nodes[source.idx].edges; i != (uint32_t) -1; i = graph.edges[i].next)
    {
        min_source = MIN(min_source, graph.edges[i].variant.start);
    } // for
    uint32_t const min_offset = min_source - graph.nodes[source.idx].row;

    graph.nodes[source.idx].row += min_offset;
    graph.nodes[source.idx].col += min_offset;
    graph.nodes[source.idx].length -= min_offset;
    graph.nodes[sink.idx].length -= graph.nodes[sink.idx].row + graph.nodes[sink.idx].length - max_sink;

    graph.supremal = (VA_Variant) {min_source, max_sink, graph.nodes[source.idx].col - shift, graph.nodes[sink.idx].col + graph.nodes[sink.idx].length - shift};
    graph.source = source.idx;

    return graph;
} // build_graph


Graph
build3(VA_Allocator const allocator,
       size_t const len_ref, char const reference[static restrict len_ref],
       size_t const len_obs, char const observed[static restrict len_obs],
       size_t const shift)
{
    VA_LCS lcs = va_edit3(allocator, len_ref, reference, len_obs, observed);

    Graph graph = {.nodes = NULL, .edges = NULL};

    if (lcs.nodes == NULL)
    {
        uint32_t const sink_idx = add_node(allocator, &graph, len_ref, len_obs, 0);
        if (len_ref == 0 && len_obs == 0)
        {
            graph.source = sink_idx;
            return graph;
        } // if

        graph.source = add_node(allocator, &graph, shift, shift, 0);
        graph.nodes[graph.source].edges = add_edge(allocator, &graph, sink_idx, GVA_NULL, (VA_Variant) {shift, shift + len_ref, 0, len_obs});
        return graph;
    } // if

    if (lcs.length == 0 || lcs.nodes == NULL)
    {
        graph.source = add_node(allocator, &graph, shift, shift, 0);
        if (len_ref == 0 && len_obs == 0)
        {
            graph.supremal = (VA_Variant) {0, 0, 0, 0};
            return graph;
        } // if
        uint32_t const sink = add_node(allocator, &graph, len_ref, len_obs, 0);
        graph.supremal = (VA_Variant) {shift, shift + len_ref, 0, len_obs};
        graph.nodes[graph.source].edges = add_edge(allocator, &graph, sink, -1, graph.supremal);
        return graph;
    } // if

    uint32_t tail = lcs.index[lcs.length - 1].tail;
    VA_LCS_Node3 sink = lcs.nodes[tail];
    if (sink.row + sink.length == len_ref + shift && sink.col + sink.length == len_obs + shift)
    {
        lcs.nodes[tail].idx = add_node(allocator, &graph, sink.row, sink.col, sink.length);
        sink.idx = lcs.nodes[tail].idx;
    } // if
    else
    {
        sink = (VA_LCS_Node3) {.row = len_ref + shift, .col = len_obs + shift, .length = 0};
        sink.idx = add_node(allocator, &graph, sink.row, sink.col, sink.length);
        tail = GVA_NULL;
    } // else

    for (uint32_t i = lcs.index[lcs.length - 1].head; i != tail; i = lcs.nodes[i].next)
    {
        eq_count += 1;
        fprintf(stderr, "(%u, %u, %u) vs (%u, %u, %u)\n", lcs.nodes[i].row, lcs.nodes[i].col, lcs.nodes[i].length, sink.row, sink.col, sink.length);

        VA_Variant const variant = {lcs.nodes[i].row + lcs.nodes[i].length, sink.row + sink.length, lcs.nodes[i].col + lcs.nodes[i].length - shift, sink.col + sink.length - shift};

        lcs.nodes[i].idx = add_node(allocator, &graph, lcs.nodes[i].row, lcs.nodes[i].col, lcs.nodes[i].length);
        graph.nodes[lcs.nodes[i].idx].edges = add_edge(allocator, &graph, sink.idx, graph.nodes[lcs.nodes[i].idx].edges, variant);
    } // for

    for (size_t i = lcs.length - 1; i >= 1; --i)
    {
        fprintf(stderr, "level %zu\n", i);
        for (uint32_t j = 0; j < lcs.length; ++j)
        {
            fprintf(stderr, "%2u: %2d  %2d\n", j, lcs.index[j].head, lcs.index[j].tail);
        } // for
        for (size_t j = 0; j < va_array_length(lcs.nodes); ++j)
        {
            fprintf(stderr, "%2zu: (%u, %u, %u)  %2d\n", j, lcs.nodes[j].row, lcs.nodes[j].col, lcs.nodes[j].length, lcs.nodes[j].next);
        } // for

        uint32_t next = -1;
        for (uint32_t j = lcs.index[i].head; j != GVA_NULL; j = next)
        {
            next = lcs.nodes[j].next;
            if (lcs.nodes[j].idx == GVA_NULL)
            {
                continue;
            } // if

            uint32_t here = -1;
            for (uint32_t k = lcs.index[i - 1].head; k != GVA_NULL; k = lcs.nodes[k].next)
            {
                if (k > j)
                {
                    fprintf(stderr, "(%u, %u, %u) skip because topo order\n", lcs.nodes[k].row, lcs.nodes[k].col, lcs.nodes[k].length);
                    continue;
                } // if

                eq_count += 1;
                fprintf(stderr, "(%u, %u, %u) vs (%u, %u, %u)\n", lcs.nodes[k].row, lcs.nodes[k].col, lcs.nodes[k].length, lcs.nodes[j].row, lcs.nodes[j].col, lcs.nodes[j].length);

                if (lcs.nodes[k].row + lcs.nodes[k].length < lcs.nodes[j].row + lcs.nodes[j].length &&
                    lcs.nodes[k].col + lcs.nodes[k].length < lcs.nodes[j].col + lcs.nodes[j].length)
                {
                    VA_Variant const variant = {lcs.nodes[k].row + lcs.nodes[k].length, lcs.nodes[j].row + lcs.nodes[j].length - 1, lcs.nodes[k].col + lcs.nodes[k].length - shift, lcs.nodes[j].col + lcs.nodes[j].length - 1 - shift};

                    if (lcs.nodes[k].incoming == i)
                    {
                        fprintf(stderr, "SPLIT at %zu\n", i);
                        uint32_t const split_idx = lcs.nodes[k].idx;
                        lcs.nodes[k].idx = add_node(allocator, &graph, lcs.nodes[k].row, lcs.nodes[k].col, lcs.nodes[k].length);
                        lcs.nodes[k].incoming = 0;

                        // lambda-edge
                        graph.nodes[lcs.nodes[k].idx].lambda = split_idx;

                        graph.nodes[split_idx].row += lcs.nodes[k].length;
                        graph.nodes[split_idx].col += lcs.nodes[k].length;
                        graph.nodes[split_idx].length -= lcs.nodes[k].length;
                        fprintf(stderr, "RELABEL NODE %u: (%u, %u, %u)\n", split_idx, graph.nodes[split_idx].row, graph.nodes[split_idx].col, graph.nodes[split_idx].length);
                    } // if
                    else if (lcs.nodes[k].idx == GVA_NULL)
                    {
                        lcs.nodes[k].idx = add_node(allocator, &graph, lcs.nodes[k].row, lcs.nodes[k].col, lcs.nodes[k].length);
                    } // if

                    graph.nodes[lcs.nodes[k].idx].edges = add_edge(allocator, &graph, lcs.nodes[j].idx, graph.nodes[lcs.nodes[k].idx].edges, variant);

                    here = k;
                } // if
            } // for

            if (lcs.nodes[j].length > 1)
            {
                lcs.nodes[j].length -= 1;
                if (here != GVA_NULL)
                {
                    lcs.nodes[j].incoming = i;
                    lcs.nodes[j].next = lcs.nodes[here].next;
                    lcs.nodes[here].next = j;
                } // if
                else
                {
                    lcs.nodes[j].next = lcs.index[i - 1].head;
                    lcs.index[i - 1].head = j;
                } // else
                fprintf(stderr, "INSERT (%u %u %u) here: %2d\n", lcs.nodes[j].row, lcs.nodes[j].col, lcs.nodes[j].length, here);
            } // if
        } // for
    } // for

    size_t head = lcs.index[0].head;
    VA_LCS_Node3 source = lcs.nodes[head];
    if (source.row == shift && source.col == shift)
    {
        head = source.next;
    } // if
    else
    {
        source = (VA_LCS_Node3) {.row = shift, .col = shift, .length = 0};
        source.idx = add_node(allocator, &graph, source.row, source.col, source.length);
    } // else

    for (uint32_t i = head; i != GVA_NULL; i = lcs.nodes[i].next)
    {
        if (lcs.nodes[i].idx == GVA_NULL)
        {
            continue;
        } // if

        eq_count += 1;
        fprintf(stderr, "(%u, %u, %u) vs (%u, %u, %u)\n", source.row, source.col, source.length, lcs.nodes[i].row, lcs.nodes[i].col, lcs.nodes[i].length);

        VA_Variant const variant = {source.row, lcs.nodes[i].row + lcs.nodes[i].length - 1, source.col - shift, lcs.nodes[i].col + lcs.nodes[i].length - 1 - shift};

        graph.nodes[source.idx].edges = add_edge(allocator, &graph, lcs.nodes[i].idx, graph.nodes[source.idx].edges, variant);
    } // for

    uint32_t min_source = GVA_NULL;
    for (uint32_t i = graph.nodes[source.idx].edges; i != (uint32_t) -1; i = graph.edges[i].next)
    {
        min_source = MIN(min_source, graph.edges[i].variant.start);
    } // for

    if (min_source != GVA_NULL)
    {
        uint32_t const min_offset = min_source - graph.nodes[source.idx].row;
        graph.nodes[source.idx].row += min_offset;
        graph.nodes[source.idx].col += min_offset;
        graph.nodes[source.idx].length -= min_offset;
    } // if
    else
    {
        graph.nodes[source.idx].length = 0;
    } // else

    graph.source = source.idx;

    va_std_allocator.alloc(va_std_allocator.context, lcs.index, lcs.length * sizeof(*lcs.index), 0);
    va_array_destroy(va_std_allocator, lcs.nodes);
    return graph;
} // build3


void
to_dot(Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    fprintf(stderr, "digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\nsi[label=\"\",shape=none,width=0]\nsi->s%u\n", graph.source);
    for (uint32_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        fprintf(stderr, "s%u[label=\"(%u, %u, %u)\"%s]\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, graph.nodes[i].edges == (uint32_t) -1 ? ",peripheries=2" : "");
        if (graph.nodes[i].lambda != (uint32_t) -1)
        {
            fprintf(stderr, "s%u->s%u[label=\"&lambda;\",style=\"dashed\"]\n", i, graph.nodes[i].lambda);
        } // if
        for (uint32_t j = graph.nodes[i].edges; j != (uint32_t) -1; j = graph.edges[j].next)
        {
            fprintf(stderr, "s%u->s%u[label=\"%u:%u/%.*s\"]\n", i, graph.edges[j].tail, graph.edges[j].variant.start, graph.edges[j].variant.end, (int) graph.edges[j].variant.obs_end - graph.edges[j].variant.obs_start, observed + graph.edges[j].variant.obs_start);
        } // for
    } // for
    fprintf(stderr, "}\n");
} // to_dot


static size_t
bfs_traversal(VA_Allocator const allocator, Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    struct
    {
        uint32_t depth;
        uint32_t next;
    }* table = allocator.alloc(allocator.context, NULL, 0, sizeof(*table) * va_array_length(graph.nodes));

    if (table == NULL)
    {
        return -1;
    }

    for (size_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        table[i].depth = 0;
        table[i].next = -1;
    } // for
    uint32_t head = graph.source;
    uint32_t tail = graph.source;

    size_t count = 0;
    //printf("digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\nsi[shape=point,width=.1]\n");
    //printf("si->s%u\n", graph.source);

    printf("    \"edges\": [\n");
    while (head != (uint32_t) -1)
    {
        //printf("pop %u\n", head);
        for (uint32_t i = head; i != (uint32_t) -1; i = graph.nodes[i].lambda)
        {
            //printf("s%u[label=\"(%u, %u, %u)\"%s]\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, graph.nodes[i].edges == (uint32_t) -1 ? ",peripheries=2" : "");
            if (i != head)
            {
                //printf("lambda %u\n", i);
            } // if
            for (uint32_t j = graph.nodes[i].edges; j != (uint32_t) -1; j = graph.edges[j].next)
            {
                if (count > 0)
                {
                    printf(",\n");
                } // if
                count += 1;
                printf("         {\"head\": \"s%u\", \"tail\": \"s%u\", \"variant\": \"%u:%u/%.*s\"}", head, graph.edges[j].tail, graph.edges[j].variant.start, graph.edges[j].variant.end, (int) graph.edges[j].variant.obs_end - graph.edges[j].variant.obs_start, observed + graph.edges[j].variant.obs_start);
                //printf("s%u->s%u[label=\"%u:%u/%.*s\"]\n", head, graph.edges[j].tail, graph.edges[j].variant.start, graph.edges[j].variant.end, (int) graph.edges[j].variant.obs_end - graph.edges[j].variant.obs_start, observed + graph.edges[j].variant.obs_start);

                if (table[graph.edges[j].tail].depth > 0)
                {
                    //printf("skip %u\n", graph.edges[j].tail);
                    continue;
                } // if

                //printf("push %u\n", graph.edges[j].tail);
                table[graph.edges[j].tail].depth = table[i].depth + 1;
                table[tail].next = graph.edges[j].tail;
                tail = graph.edges[j].tail;
            } // for
        } // for
        head = table[head].next;

        /*
        printf("  # \tdepth\tnext\n");
        for (size_t i = 0; i < va_array_length(graph.nodes); ++i)
        {
            printf("%3zu:\t%5u\t%4d\n", i, table[i].depth, table[i].next);
        } // for
        printf("head: %d\ntail: %d\n", head, tail);
        */
    } // while

    printf("\n    ]\n");

    table = allocator.alloc(allocator.context, table, sizeof(*table) * va_array_length(graph.nodes), 0);

    return count;
} // bfs_traversal


static void
lambda_edges(Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    printf("    \"edges\": [\n");
    size_t count = 0;
    for (uint32_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        if (graph.nodes[i].lambda != (uint32_t) -1)
        {
            if (count > 0)
            {
                printf(",\n");
            } // if
            count += 1;
            printf("         {\"head\": \"s%u\", \"tail\": \"s%u\", \"variant\": \"lambda\"}", i, graph.nodes[i].lambda);
        } // if
        for (uint32_t j = graph.nodes[i].edges; j != (uint32_t) -1; j = graph.edges[j].next)
        {
            if (count > 0)
            {
                printf(",\n");
            } // if
            count += 1;
            printf("         {\"head\": \"s%u\", \"tail\": \"s%u\", \"variant\": \"%u:%u/%.*s\"}", i, graph.edges[j].tail, graph.edges[j].variant.start, graph.edges[j].variant.end, (int) graph.edges[j].variant.obs_end - graph.edges[j].variant.obs_start, observed + graph.edges[j].variant.obs_start);
        } // for
    } // for
    printf("\n    ]\n");
} // lambda_edges


void
to_json(Graph const graph, size_t const len_obs, char const observed[static len_obs], bool const lambda)
{
    printf("{\n    \"source\": \"s%u\",\n    \"nodes\": {\n", graph.source);
    for (uint32_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        printf("        \"s%u\": {\"row\": %u, \"col\": %u, \"length\": %u}%s\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, i < va_array_length(graph.nodes) - 1 ? "," : "");
    } // for
    printf("    },\n");
    if (lambda) {
        lambda_edges(graph, len_obs, observed);
    } else {
        bfs_traversal(va_std_allocator, graph, len_obs, observed);
    }
    //printf("    \"supremal\": \"%u:%u/%.*s\",\n", graph.supremal.start, graph.supremal.end, (int) graph.supremal.obs_end - graph.supremal.obs_start, observed + graph.supremal.obs_start);
    //printf("    \"local_supremal\": [\n");
    //local_supremal(va_std_allocator, graph, len_obs, observed, false);
    //printf("\n    ],\n");
    //printf("    \"canonical\": [\n");
    //canonical(va_std_allocator, graph, len_obs, observed, false);
    //printf("\n    ]\n");
    printf("}");
} // to_json


size_t
edges(uint32_t const head_row, uint32_t const head_col, uint32_t const head_length,
      uint32_t const tail_row, uint32_t const tail_col, uint32_t const tail_length,
      bool const is_source,
      size_t const len_obs, char const observed[static len_obs],
      VA_Variant const edge)
{
    printf(
        "(%u, %u, %u) -> (%u, %u, %u) source: %d\n",
        head_row, head_col, head_length,
        tail_row, tail_col, tail_length,
        is_source
    );
    //printf(VAR_FMT "\n", print_variant(edge, observed));
    printf("%u %u %u %u\n", edge.start, edge.end, edge.obs_start, edge.obs_end);

    ptrdiff_t const row = (ptrdiff_t) head_row - is_source;
    ptrdiff_t const col = (ptrdiff_t) head_col - is_source;
    uint32_t const length = head_length + is_source;

    printf("%zd, %zd, %u\n", row, col, length);

    ptrdiff_t const offset = MIN((ptrdiff_t) tail_row - row, (ptrdiff_t) tail_col - col) - 1;
    uint32_t const head_start = offset > 0 ? MIN(offset, length - 1) : 0;
    uint32_t const tail_start = offset < 0 ? MIN(-offset, tail_length - 1) : 0;

    uint32_t const count = MIN(length - head_start, tail_length - tail_start);

    printf("    %zd: (%u, %u) :: %u\n", offset, head_start, tail_start, count);

    if (offset < 0 && tail_length <= -offset)
    {
        printf("    X\n");
        return 0;
    } // if

    bool found = false;
    for (uint32_t j = 0; j < count; ++j)
    {
        VA_Variant const variant =
        {
            row + head_start + j + 1,
            tail_row + tail_start + j,
            col + head_start + j + 1,
            tail_col + tail_start + j
        };

        if (edge.start == variant.start && edge.end == variant.end &&
            edge.obs_start == variant.obs_start && edge.obs_end == variant.obs_end)
        {
            found = true;
        } // if

        printf("    (%u, %u) -> (%u, %u) :: " VAR_FMT "\n", //"%u %u %u %u\n",
            head_row + head_start + j, head_col + head_start + j,
            tail_row + tail_start + j, tail_col + tail_start + j,
            print_variant(variant, observed) //,
            //variant.start, variant.end, variant.obs_start, variant.obs_end
        );
    } // for

    if (!found)
    {
        return -1;
    } // if
    return count;
} // edges
