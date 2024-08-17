#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint32_t
#include <stdio.h>      // stderr, fprintf, printf
#include <stdlib.h>     // EXIT_*
#include <string.h>     // strlen

#include "../include/alloc.h"           // VA_Allocator
#include "../include/array.h"           // va_array_*
#include "../include/edit.h"            // va_edit
#include "../include/std_alloc.h"       // va_std_allocator
#include "../include/variant.h"         // VA_Variant


#define MAX(lhs, rhs) (((lhs) > (rhs)) ? (lhs) : (rhs))
#define MIN(lhs, rhs) (((lhs) < (rhs)) ? (lhs) : (rhs))

#define print_variant(variant, observed) variant.start, variant.end, (int) variant.obs_end - variant.obs_start, observed + variant.obs_start
#define VAR_FMT "%u:%u/%.*s"


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
    VA_Variant supremal;
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


static inline void
destroy(VA_Allocator const allocator, Graph* const graph)
{
    graph->nodes = va_array_destroy(allocator, graph->nodes);
    graph->edges = va_array_destroy(allocator, graph->edges);
} // destroy


static Graph
build_graph(VA_Allocator const allocator,
            size_t const len_ref, size_t const len_obs,
            size_t const len_lcs, VA_LCS_Node* lcs_nodes[static len_lcs],
            size_t const shift)
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
    size_t here = 0;
    VA_LCS_Node* const heads = lcs_nodes[len_lcs - 1];
    for (size_t i = 0; i < va_array_length(heads); ++i)
    {
        if (heads[i].row + heads[i].length < sink.row + sink.length && heads[i].col + heads[i].length < sink.col + sink.length)
        {
            VA_Variant const variant = {heads[i].row + heads[i].length, sink.row + sink.length - 1, heads[i].col + heads[i].length - shift, sink.col + sink.length - 1 - shift};
            max_sink = MAX(max_sink, sink.row + sink.length - 1);
            heads[i].idx = add_node(allocator, &graph, heads[i].row, heads[i].col, heads[i].length);
            graph.nodes[heads[i].idx].edges = add_edge(allocator, &graph, sink.idx, graph.nodes[heads[i].idx].edges, variant);
            here = i + 1;
        } // if
    } // for

    if (sink.length > 1)
    {
        sink.length -= 1;
        sink.incoming = len_lcs - 1;
        va_array_insert(allocator, lcs_nodes[len_lcs - 1], sink, here);
    } // if

    for (size_t i = len_lcs - 1; i >= 1; --i)
    {
        for (size_t j = 0; j < va_array_length(lcs_nodes[i]); ++j)
        {
            if (lcs_nodes[i][j].idx == (uint32_t) -1)
            {
                continue;
            } // if

            here = 0;
            VA_LCS_Node const tail = lcs_nodes[i][j];
            VA_LCS_Node* const heads = lcs_nodes[i - 1];
            for (size_t k = 0; k < va_array_length(heads); ++k)
            {
                if (heads[k].row + heads[k].length < tail.row + tail.length && heads[k].col + heads[k].length < tail.col + tail.length)
                {
                    VA_Variant const variant = {heads[k].row + heads[k].length, tail.row + tail.length - 1, heads[k].col + heads[k].length - shift, tail.col + tail.length - 1 - shift};
                    max_sink = MAX(max_sink, tail.row + tail.length - 1);

                    if (heads[k].incoming == i)
                    {
                        size_t const split_idx = heads[k].idx;
                        heads[k].idx = add_node(allocator, &graph, heads[k].row, heads[k].col, graph.nodes[split_idx].length);
                        heads[k].incoming = 0;

                        // lambda-edge
                        graph.nodes[heads[k].idx].lambda = split_idx;

                        graph.nodes[heads[k].idx].edges = add_edge(allocator, &graph, tail.idx, graph.nodes[heads[k].idx].edges, variant);

                        graph.nodes[split_idx].row += heads[k].length;
                        graph.nodes[split_idx].col += heads[k].length;
                        graph.nodes[split_idx].length -= heads[k].length;
                    } // if
                    else
                    {
                        if (heads[k].idx == (uint32_t) -1)
                        {
                            heads[k].idx = add_node(allocator, &graph, heads[k].row, heads[k].col, heads[k].length);
                        } // if
                        graph.nodes[heads[k].idx].edges = add_edge(allocator, &graph, tail.idx, graph.nodes[heads[k].idx].edges, variant);
                    } // else
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
            max_sink = MAX(max_sink, lcs_nodes[0][i].row + lcs_nodes[0][i].length - 1);
            graph.nodes[source.idx].edges = add_edge(allocator, &graph, lcs_nodes[0][i].idx, graph.nodes[source.idx].edges, variant);
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


static void
to_dot(Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    fprintf(stderr, "digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\nsi[shape=point,width=.1]\nsi->s%u\n", graph.source);
    for (uint32_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        fprintf(stderr, "s%u[label=\"(%u, %u, %u)\"%s]\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, graph.nodes[i].edges == (uint32_t) -1 ? ",peripheries=2" : "");
        if (graph.nodes[i].lambda != (uint32_t) -1)
        {
            fprintf(stderr, "s%u->s%u[label=\"&#955;\",style=\"dashed\"]\n", i, graph.nodes[i].lambda);
        } // if
        for (uint32_t j = graph.nodes[i].edges; j != (uint32_t) -1; j = graph.edges[j].next)
        {
            fprintf(stderr, "s%u->s%u[label=\"%u:%u/%.*s\"]\n", i, graph.edges[j].tail, graph.edges[j].variant.start, graph.edges[j].variant.end, (int) graph.edges[j].variant.obs_end - graph.edges[j].variant.obs_start, observed + graph.edges[j].variant.obs_start);
        } // for
    } // for
    fprintf(stderr, "}\n");
} // to_dot


typedef struct
{
    uint32_t post;
    uint32_t rank;
    uint32_t edge;
    uint32_t start;
    uint32_t end;
    uint32_t next;
} Post_Table;


static uint32_t
intersection(uint32_t lhs, uint32_t rhs, size_t const length, Post_Table const visited[static length])
{
    if (lhs == (uint32_t) -1)
    {
        return rhs;
    } // if
    while (lhs != rhs)
    {
        while (visited[lhs].rank > visited[rhs].rank)
        {
            lhs = visited[lhs].post;
        } // while
        while (visited[rhs].rank > visited[lhs].rank)
        {
            rhs = visited[rhs].post;
        } // while
    } // while
    return lhs;
} // intersection


static void
local_supremal(VA_Allocator const allocator, Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    uint32_t const length = va_array_length(graph.nodes);
    Post_Table* visited = allocator.alloc(allocator.context, NULL, 0, sizeof(*visited) * length);
    if (visited == NULL)
    {
        return;
    } // if

    for (uint32_t i = 0; i < length; ++i)
    {
        visited[i].post = -1;
        visited[i].start = 0;
        visited[i].end = -1;
        visited[i].next = -1;
    } // for

    uint32_t rank = 0;
    uint32_t head = graph.source;
    visited[head].edge = graph.nodes[head].edges;

    while (head != (uint32_t) -1)
    {
        uint32_t edge = visited[head].edge;
        if (edge == (uint32_t) -1)
        {
            uint32_t const lambda = graph.nodes[head].lambda;
            if (lambda != (uint32_t) -1)
            {
                if (visited[lambda].next == (uint32_t) -1)
                {
                    visited[lambda].edge = graph.nodes[lambda].edges;
                    visited[lambda].next = head;
                    head = lambda;
                    continue;
                } // if
                visited[head].post = intersection(visited[head].post, lambda, length, visited);
            } // if

            visited[head].rank = rank;
            rank += 1;
            head = visited[head].next;
            continue;
        } // if

        if (visited[graph.edges[edge].tail].next != (uint32_t) -1)
        {
            visited[head].post = intersection(visited[head].post, graph.edges[edge].tail, length, visited);
            visited[head].edge = graph.edges[edge].next;
        } // if
        while (edge != (uint32_t) -1)
        {
            uint32_t const tail = graph.edges[edge].tail;
            visited[head].end = MIN(visited[head].end, graph.edges[edge].variant.start);
            visited[tail].start = MAX(visited[tail].start, graph.edges[edge].variant.end);
            if (visited[tail].next == (uint32_t) -1)
            {
                visited[tail].edge = graph.nodes[tail].edges;
                visited[tail].next = head;
                head = tail;
                break;
            } // if
            visited[head].post = intersection(visited[head].post, tail, length, visited);
            edge = graph.edges[edge].next;
        } // while
    } // while

    /*
    printf(" # \tpost\trank\tstart\tend\n");
    for (uint32_t i = 0; i < length; ++i)
    {
        printf("%2u:\t%4d\t%4u\t%5u\t%3d\n", i, visited[i].post, visited[i].rank, visited[i].start, visited[i].end);
    } // for
    */

    head = graph.source;
    uint32_t shift = 0;
    for (uint32_t tail = visited[head].post; tail != (uint32_t) -1; head = tail, tail = visited[head].post)
    {
        uint32_t const start = visited[head].end;
        uint32_t const end = visited[tail].start;
        VA_Variant const variant = {start, end, graph.nodes[head].col + start - graph.nodes[head].row - shift, graph.nodes[tail].col + end - graph.nodes[tail].row - shift};
        if (head != graph.source)
        {
            printf(",\n");
        } // if
        printf("        \"" VAR_FMT "\"", print_variant(variant, observed));
    } // for

    visited = allocator.alloc(allocator.context, visited, sizeof(*visited) * length, 0);
} // local_supremal


typedef struct
{
    uint32_t lca;
    uint32_t rank;
    uint32_t depth;
    uint32_t start;
    uint32_t end;
    uint32_t prev;
    uint32_t next;
} LCA_Table;


static inline uint32_t
lca(uint32_t* const start, uint32_t lhs, uint32_t rhs, size_t const length, LCA_Table const visited[static length])
{
    uint32_t lhs_start = *start;
    uint32_t rhs_start = *start;
    while (lhs != rhs)
    {
        fprintf(stderr, "        lca %u (%d), %u (%d)\n", lhs, lhs_start, rhs, rhs_start);

        while (visited[lhs].rank > visited[rhs].rank)
        {
            fprintf(stderr, "        left\n");
            lhs_start = visited[lhs].start;
            lhs = visited[lhs].lca;
        } // if
        while (visited[rhs].rank > visited[lhs].rank)
        {
            fprintf(stderr, "        right\n");
            rhs_start = visited[rhs].start;
            rhs = visited[rhs].lca;
        } // if
    } // while
    *start = MIN(lhs_start, rhs_start);
    return lhs;
} // lca


static void
canonical(VA_Allocator const allocator, Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    uint32_t const length = va_array_length(graph.nodes);
    LCA_Table* visited = allocator.alloc(allocator.context, NULL, 0, sizeof(*visited) * length);
    if (visited == NULL)
    {
        return;
    } // if

    for (uint32_t i = 0; i < length; ++i)
    {
        visited[i].depth = -1;
    } // for

    uint32_t sink = -1;
    visited[graph.source] = (LCA_Table) {-1, 0, 0, -1, -1, -1, -1};
    uint32_t rank = 1;
    // main loop over a queue
    for (uint32_t head = graph.source, tail = graph.source; head != (uint32_t) -1; head = visited[head].next)
    {
        fprintf(stderr, "pop %u @ %u\n", head, visited[head].depth);

        if (graph.nodes[head].edges == (uint32_t) -1)
        {
            fprintf(stderr, "    sink\n");
            sink = head;
            continue;
        } // if

        uint32_t const lambda = graph.nodes[head].lambda;
        if (lambda != (uint32_t) -1)
        {
            if (visited[lambda].depth == (uint32_t) -1)
            {
                // add lambda in stack order
                fprintf(stderr, "    push lambda %u @ %u\n", lambda, visited[head].depth);
                visited[lambda] = (LCA_Table) {head, rank, visited[head].depth, -1, -1, head, visited[head].next};
                rank += 1;
                visited[head].next = lambda;
            } // if
            else if (visited[lambda].depth == visited[head].depth)
            {
                fprintf(stderr, "    update lambda %u @ %u\n", lambda, visited[head].depth);
                fprintf(stderr, "        lca(%u, %u)\n", visited[lambda].lca, head);
                visited[lambda].lca = lca(&visited[lambda].start, visited[lambda].lca, head, length, visited);

                // FIXME: what is correct here?
                // AGGCCG CGCAGCCTC        ignores the common end position
                // TAATTGGTAG CTATTCAG     needs the common end position

                //visited[lambda].end = graph.nodes[head].row;  // this simulates the empty variant on a lambda edge

                fprintf(stderr, "        = %u (%u:%u)\n", visited[lambda].lca, visited[lambda].start, visited[lambda].end);
            } // if
            else
            {
                fprintf(stderr, "    improve lambda %u @ %u\n", lambda, visited[head].depth);

                // FIXME: what about the original link to lambda?
                //        lambda is already in the queue and should be *moved* instead of added again (in stack order)
                //        if moving is difficult we can possibly skip nodes that are already processed (in the main loop)
                // possibly related:
                // INFINITE LOOP: CAGCGAGT CGTGTAAGGTGTACTGAAA
                // FIXED by double linked list

                visited[visited[lambda].prev].next = visited[lambda].next;
                visited[lambda] = (LCA_Table) {head, visited[lambda].rank, visited[head].depth, -1, -1, head, visited[head].next};
                visited[head].next = lambda;
            } // else
        } // if

        for (uint32_t i = graph.nodes[head].edges; i != (uint32_t) -1; i = graph.edges[i].next)
        {
            uint32_t const edge_tail = graph.edges[i].tail;
            if (visited[edge_tail].depth == (uint32_t) -1)
            {
                // add regular successors in queue order
                fprintf(stderr, "    push %u @ %u\n", edge_tail, visited[head].depth + 1);
                visited[edge_tail] = (LCA_Table) {head, rank, visited[head].depth + 1, graph.edges[i].variant.start, graph.edges[i].variant.end, tail, -1};
                rank += 1;
                visited[tail].next = edge_tail;
                tail = edge_tail;
            } // if
            else if (visited[edge_tail].depth - 1 == visited[head].depth)
            {
                fprintf(stderr, "    update %u @ %u\n", edge_tail, visited[edge_tail].depth);
                fprintf(stderr, "        lca(%u, %u)\n", visited[edge_tail].lca, head);
                visited[edge_tail].start = MIN(visited[edge_tail].start, graph.edges[i].variant.start);
                visited[edge_tail].lca = lca(&visited[edge_tail].start, visited[edge_tail].lca, head, length, visited);
                visited[edge_tail].end = MAX(visited[edge_tail].end, graph.edges[i].variant.end);
                fprintf(stderr, "        = %u (%u:%u)\n", visited[edge_tail].lca, visited[edge_tail].start, visited[edge_tail].end);
            } // if
            else
            {
                fprintf(stderr, "    skip %u @ %u\n", edge_tail, visited[edge_tail].depth);
            } // else
        } // for
    } // for

    fprintf(stderr, " # \tlca\trank\tdepth\tstart\tend\tprev\tnext\n");
    for (uint32_t i = 0; i < length; ++i)
    {
        fprintf(stderr, "%2u:\t%3d\t%4u\t%5u\t%5d\t%3d\t%4d\t%4d\n", i, visited[i].lca, visited[i].rank, visited[i].depth, visited[i].start, visited[i].end, visited[i].prev, visited[i].next);
    } // for

    // the canonical variant is given in the reverse order (for now): the Python checker needs to reverse this list
    for (uint32_t tail = sink; tail != (uint32_t) -1; tail = visited[tail].lca)
    {
        uint32_t const head = visited[tail].lca;
        if (head != (uint32_t) - 1 && visited[tail].start != (uint32_t) -1)
        {
            // we need to skip lambda edges: start && end == -1
            uint32_t const start_offset = visited[tail].start - graph.nodes[head].row;
            uint32_t const end_offset = visited[tail].end - graph.nodes[tail].row;
            VA_Variant const variant =  {visited[tail].start, visited[tail].end, graph.nodes[head].col + start_offset, graph.nodes[tail].col + end_offset};

            if (tail != sink)
            {
                printf(",\n");
            } // if
            printf("        \"" VAR_FMT "\"", print_variant(variant, observed));
        } // if
    } // while

    visited = allocator.alloc(allocator.context, visited, sizeof(*visited) * length, 0);
} // canonical


static size_t
bfs_traversal(VA_Allocator const allocator, Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    struct
    {
        uint32_t depth;
        uint32_t next;
    }* table = allocator.alloc(allocator.context, NULL, 0, sizeof(*table) * va_array_length(graph.nodes));

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

    printf("\n    ],\n");

    table = allocator.alloc(allocator.context, table, sizeof(*table) * va_array_length(graph.nodes), 0);

    return count;
} // bfs_traversal


static void
to_json(Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    printf("{\n    \"source\": \"s%u\",\n    \"nodes\": {\n", graph.source);
    for (uint32_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        printf("        \"s%u\": {\"row\": %u, \"col\": %u, \"length\": %u}%s\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, i < va_array_length(graph.nodes) - 1 ? "," : "");
    } // for
    printf("    },\n");
    bfs_traversal(va_std_allocator, graph, len_obs, observed);
    printf("    \"supremal\": \"%u:%u/%.*s\",\n", graph.supremal.start, graph.supremal.end, (int) graph.supremal.obs_end - graph.supremal.obs_start, observed + graph.supremal.obs_start);
    printf("    \"local_supremal\": [\n");
    local_supremal(va_std_allocator, graph, len_obs, observed);
    printf("\n    ],\n");
    printf("    \"canonical\": [\n");
    canonical(va_std_allocator, graph, len_obs, observed);
    printf("    ]\n");
    printf("}");
} // to_json


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
    */

    Graph graph = build_graph(va_std_allocator, len_ref, len_obs, len_lcs, lcs_nodes, 0);

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

    to_dot(graph, len_obs, observed);

    //local_supremal(va_std_allocator, graph, len_obs, observed);

    //canonical(va_std_allocator, graph, len_obs, observed);

    to_json(graph, len_obs, observed);

    destroy(va_std_allocator, &graph);

    return EXIT_SUCCESS;
} // main
