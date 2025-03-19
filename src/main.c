#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // intmax_t, uint32_t
#include <stdio.h>      // stderr, fprintf, printf
#include <stdlib.h>     // EXIT_*, atoi, rand, srand
#include <string.h>     // strlen, strncmp

#include "../include/alloc.h"           // VA_Allocator
#include "../include/array.h"           // va_array_*
#include "../include/edit.h"            // va_edit
#include "../include/std_alloc.h"       // va_std_allocator
#include "../include/variant.h"         // VA_Variant


#define MAX(lhs, rhs) (((lhs) > (rhs)) ? (lhs) : (rhs))
#define MIN(lhs, rhs) (((lhs) < (rhs)) ? (lhs) : (rhs))

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


typedef struct
{
    uint32_t tail;
    uint32_t next;
    VA_Variant variant;
} Edge;


typedef struct
{
    uint32_t tail;
    uint32_t next;
} Edge2;


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


typedef struct
{
    Node* nodes;
    Edge2* edges;
    uint32_t source;
} Graph2;


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
            //printf("MAKE EDGE (%u %u %u) -> (%u %u %u)\n", heads[i].row, heads[i].col, heads[i].length, sink.row, sink.col, sink.length);

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

        //printf("INSERT (%u %u %u) here: %zu\n", sink.row, sink.col, sink.length, here);
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
                    max_sink = MAX(max_sink, tail.row + tail.length - 1);

                    //printf("MAKE EDGE (%u %u %u) -> (%u %u %u)\n", heads[k].row, heads[k].col, heads[k].length, tail.row, tail.col, tail.length);

                    VA_Variant const variant = {heads[k].row + heads[k].length, tail.row + tail.length - 1, heads[k].col + heads[k].length - shift, tail.col + tail.length - 1 - shift};

                    if (heads[k].incoming == i)
                    {
                        size_t const split_idx = heads[k].idx;
                        heads[k].idx = add_node(allocator, &graph, heads[k].row, heads[k].col, heads[k].length);
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

                //printf("INSERT (%u %u %u) here: %zu\n", lcs_nodes[i][j].row, lcs_nodes[i][j].col, lcs_nodes[i][j].length, here);

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
            max_sink = MAX(max_sink, lcs_nodes[0][i].row + lcs_nodes[0][i].length - 1);

            //printf("MAKE EDGE (%u %u %u) -> (%u %u %u)\n", source.row, source.col, source.length, lcs_nodes[0][i].row, lcs_nodes[0][i].col, lcs_nodes[0][i].length);

            VA_Variant const variant = {source.row, lcs_nodes[0][i].row + lcs_nodes[0][i].length - 1, source.col - shift, lcs_nodes[0][i].col + lcs_nodes[0][i].length - 1 - shift};
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
print_lca_table(size_t const length, LCA_Table const visited[static length], uint32_t const head, uint32_t const tail)
{
    fprintf(stderr, " # \tlca\trank\tdepth\tstart\tend\tprev\tnext\n");
    for (uint32_t i = 0; i < length; ++i)
    {
        fprintf(stderr, "%2u:\t%3d\t%4d\t%5d\t%5d\t%3d\t%4d\t%4d\n", i, visited[i].lca, visited[i].rank, visited[i].depth, visited[i].start, visited[i].end, visited[i].prev, visited[i].next);
    } // for
    fprintf(stderr, "head: %d\ntail: %d\n", head, tail);
} // print_lca_table


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

        visited[i].lca = -1;
        visited[i].rank = -1;
        visited[i].start = -1;
        visited[i].end = -1;
        visited[i].prev = -1;
        visited[i].next = -1;

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

                visited[visited[head].next].prev = lambda;
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
                // ACAGGA CAAGGCG
                // TAATTGGTAG CTATTCAG     needs the common end position
                // TTAAG TAACA

                //visited[lambda].end = graph.nodes[head].row;  // this simulates the empty variant on a lambda edge

                fprintf(stderr, "        = %u (%u:%u)\n", visited[lambda].lca, visited[lambda].start, visited[lambda].end);
            } // if
            else if (visited[head].next == lambda)
            {
                fprintf(stderr, "    improve lambda in place %u @ %u\n", lambda, visited[head].depth);
                visited[lambda] = (LCA_Table) {head, visited[lambda].rank, visited[head].depth, -1, -1, visited[lambda].prev, visited[lambda].next};
            } // if
            else
            {
                fprintf(stderr, "    improve lambda %u @ %u\n", lambda, visited[head].depth);

                if (lambda == tail)
                {
                    tail = visited[lambda].prev;
                } // if
                else
                {
                    visited[visited[lambda].next].prev = visited[lambda].prev;
                } // else
                visited[visited[head].next].prev = lambda;
                visited[visited[lambda].prev].next = visited[lambda].next;
                visited[lambda] = (LCA_Table) {head, visited[lambda].rank, visited[head].depth, -1, -1, head, visited[head].next};
                visited[head].next = lambda;
            } // else

            //print_lca_table(length, visited, head, tail);

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

            //print_lca_table(length, visited, head, tail);
        } // for
    } // for

    print_lca_table(length, visited, -1, -1);

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
    printf("\n    ],\n");
} // lambda_edges

/*
static void
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
    printf("    \"supremal\": \"%u:%u/%.*s\",\n", graph.supremal.start, graph.supremal.end, (int) graph.supremal.obs_end - graph.supremal.obs_start, observed + graph.supremal.obs_start);
    printf("    \"local_supremal\": [\n");
    local_supremal(va_std_allocator, graph, len_obs, observed);
    printf("\n    ],\n");
    printf("    \"canonical\": [\n");
    canonical(va_std_allocator, graph, len_obs, observed);
    printf("\n    ]\n");
    printf("}");
} // to_json
*/


static size_t
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

    Graph graph = build_graph(va_std_allocator, len_ref, len_obs, len_lcs, lcs_nodes, 0);

    if (debug)
    {
        to_dot(graph, len_obs, observed);
    } // if

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


static uint32_t
edges2(VA_LCS_Node const head, VA_LCS_Node const tail,
       bool const is_source, bool const is_sink,
       size_t const len_obs, char const observed[static len_obs], uint32_t* const start, uint32_t* const end)
{
//    printf("%u %u %u -> %u %u %u %s %s\n", head.row, head.col, head.length, tail.row, tail.col, tail.length, is_source ? "SOURCE" : "", is_sink ? "SINK" : "");

    intmax_t const row = (intmax_t) head.row - is_source;
    intmax_t const col = (intmax_t) head.col - is_source;
    uint32_t const head_length = head.length + is_source;
    uint32_t const tail_length = tail.length + is_sink;

    intmax_t const offset = MIN((intmax_t) tail.row - row, (intmax_t) tail.col - col) - 1;

    uint32_t const head_offset = offset > 0 ? MIN(head_length, offset + 1) : 1;
    uint32_t const tail_offset = offset < 0 ? MIN(tail_length, -offset) : 0;

    //printf("    %ld: (%u, %u)\n", offset, head_offset, tail_offset);

    if (head_offset > head_length || (tail_length > 0 && tail_offset >= tail_length))
    {
        printf("NO EDGE\n");
        return 0;
    } // if

    uint32_t const count = MIN(head_length - head_offset, tail_length - tail_offset - 1) + 1;

    *start = row + head_offset;
    *end = tail.row + tail_offset;
    VA_Variant const variant = {*start, *end, col + head_offset, tail.col + tail_offset};

    printf(VAR_FMT " x %u\n", print_variant(variant, observed), count);
    return count;
} // edges2


static size_t
bfs_traversal2(Graph2 const graph, size_t const len_obs, char const observed[static len_obs])
{
    struct
    {
        uint32_t depth;
        uint32_t next;
    }* table = va_std_allocator.alloc(va_std_allocator.context, NULL, 0, sizeof(*table) * va_array_length(graph.nodes));

    for (size_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        table[i].depth = 0;
        table[i].next = GVA_NULL;
    } // for
    uint32_t head = graph.source;
    uint32_t tail = graph.source;

    size_t count = 0;
    //printf("digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\nsi[shape=point,width=.1]\n");
    //printf("si->s%u\n", graph.source);

    printf("    \"edges\": [\n");
    while (head != GVA_NULL)
    {
        //printf("pop %u\n", head);
        for (uint32_t i = head; i != GVA_NULL; i = graph.nodes[i].lambda)
        {
            //printf("s%u[label=\"(%u, %u, %u)\"%s]\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, graph.nodes[i].edges == (uint32_t) -1 ? ",peripheries=2" : "");
            if (i != head)
            {
                //printf("lambda %u\n", i);
            } // if
            for (uint32_t j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
            {
                if (count > 0)
                {
                    printf(",\n");
                } // if
                count += 1;
                uint32_t const count2 = edges2(((VA_LCS_Node) {graph.nodes[head].row, graph.nodes[head].col, graph.nodes[head].length, 0, 0}),
                                               ((VA_LCS_Node) {graph.nodes[graph.edges[j].tail].row, graph.nodes[graph.edges[j].tail].col, graph.nodes[graph.edges[j].tail].length, 0, 0}),
                                               graph.nodes[head].row == graph.nodes[graph.source].row && graph.nodes[head].col == graph.nodes[graph.source].col, graph.nodes[graph.edges[j].tail].edges == GVA_NULL, len_obs, observed, &(uint32_t) {0}, &(uint32_t) {0});

                for (uint32_t k = 0; k < count2; ++k)
                {
                    printf("         {\"head\": \"s%u\", \"tail\": \"s%u\", \"variant\": \"%u\"}\n", head, graph.edges[j].tail, count2);
                } // for
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

    table = va_std_allocator.alloc(va_std_allocator.context, table, sizeof(*table) * va_array_length(graph.nodes), 0);

    return count;
} // bfs_traversal2


static Graph2
build(size_t const len_ref, char const reference[static len_ref],
      size_t const len_obs, char const observed[static len_obs],
      size_t const shift)
{
    VA_LCS_Node** lcs_nodes = NULL;
    size_t const len_lcs = va_edit(va_std_allocator, len_ref, reference, len_obs, observed, &lcs_nodes);

    Graph2 graph = {NULL, NULL, GVA_NULL};

    if (len_lcs == 0 || lcs_nodes == NULL)
    {
        fprintf(stderr, "Trivial graph\n");
        return graph;
    } // if

    VA_LCS_Node* last = &lcs_nodes[len_lcs - 1][va_array_length(lcs_nodes[len_lcs - 1]) - 1];
    if (last->row + last->length != len_ref + shift || last->col + last->length != len_obs + shift)
    {
        va_array_append(va_std_allocator, lcs_nodes[len_lcs - 1],
                        ((VA_LCS_Node) {.row = len_ref + shift, .col = len_obs + shift, .length = 0}));
        last = &lcs_nodes[len_lcs - 1][va_array_length(lcs_nodes[len_lcs - 1]) - 1];
    } // if
    va_array_append(va_std_allocator, graph.nodes, ((Node) {last->row, last->col, last->length, GVA_NULL, GVA_NULL}));
    last->idx = va_array_length(graph.nodes) - 1;
    last->incoming = -1;

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

    VA_LCS_Node source = {shift, shift, 0, -1, GVA_NULL};
    bool found_source = false;
    bool is_sink = true;
    for (size_t t_i = 0; t_i < len_lcs; ++t_i)
    {
        size_t const t_len = va_array_length(lcs_nodes[len_lcs - t_i - 1]);
        for (size_t t_j = 0; t_j < t_len; ++t_j)
        {
            VA_LCS_Node* const tail = &lcs_nodes[len_lcs - t_i - 1][t_len - t_j - 1];
            if (tail->idx == GVA_NULL)
            {
                continue;
            } // if
            printf("%u %u %u %u\n", tail->row, tail->col, tail->length, tail->incoming);

            uint32_t split_idx = tail->idx;
            for (size_t h_i = 0; h_i < MIN(len_lcs - t_i, tail->length + 1); ++h_i)
            {
                size_t const h_len = h_i == 0 ? t_len - t_j - 1: va_array_length(lcs_nodes[len_lcs - t_i - h_i - 1]);
                //printf("    level: %zu with %zu nodes\n", len_lcs - t_i - h_i - 1, h_len);
                for (size_t h_j = 0; h_j < h_len; ++h_j)
                {
                    VA_LCS_Node* head = &lcs_nodes[len_lcs - t_i - h_i - 1][h_j];

                    printf("    %u %u %u: ", head->row, head->col, head->length);

                    if (head->row + head->length + h_i > tail->row + tail->length ||
                        head->col + head->length + h_i > tail->col + tail->length)
                    {
                        if (tail->row + tail->length > head->row + head->length + h_i ||
                            tail->col + tail->length > head->col + head->length + h_i)
                        {
                            printf("--\n");  // FIXME: skip this remainder of this (head) level
                        } // if
                        else
                        {
                            printf("CONVERSE  ");
                            if (head->idx != GVA_NULL && edges2(*tail, *head, tail->row == shift && tail->col == shift,
                                                                false, len_obs, observed, &(uint32_t) {0}, &(uint32_t) {0}))
                            {
                                va_array_append(va_std_allocator, graph.edges, ((Edge2) {head->idx, graph.nodes[tail->idx].edges}));
                                graph.nodes[tail->idx].edges = va_array_length(graph.edges) - 1;
                            }
                        } // else
                    } // if
                    else
                    {
                        uint32_t start = -1;
                        uint32_t end = -1;
                        uint32_t const count = edges2(*head, *tail, head->row == shift && head->col == shift, is_sink, len_obs, observed, &start, &end);
                        if (count > 0)
                        {
                            if (head->idx == GVA_NULL)
                            {
                                va_array_append(va_std_allocator, graph.nodes, ((Node) {head->row, head->col, head->length, GVA_NULL, GVA_NULL}));
                                head->idx = va_array_length(graph.nodes) - 1;
                                if (head->row == shift && head->col == shift)
                                {
                                    source = *head;
                                    graph.source = source.idx;
                                    found_source = true;
                                } // if
                            } // if
                            head->incoming = MIN(head->incoming, start);
                            if (end + count > tail->incoming)
                            {
                                uint32_t const split_len = tail->incoming - tail->row;
                                va_array_append(va_std_allocator, graph.nodes,
                                                ((Node) {tail->row, tail->col, split_len, GVA_NULL, tail->idx}));
                                split_idx = va_array_length(graph.nodes) - 1;
                                graph.nodes[tail->idx].row += split_len;
                                graph.nodes[tail->idx].col += split_len;
                                graph.nodes[tail->idx].length -= split_len;
                                tail->incoming = end + count;
                                printf("SPLIT %u %u %u\n", tail->row, tail->col, tail->length);
                            } // if
                            va_array_append(va_std_allocator, graph.edges, ((Edge2) {tail->idx, graph.nodes[head->idx].edges}));
                            graph.nodes[head->idx].edges = va_array_length(graph.edges) - 1;
                        } // if
                        else
                        {
                           // FIXME: skip this remainder of this (head) level
                        } // else
                    } // else
                } // for h_j
                is_sink = false;
                tail->idx = split_idx;
            } // for h_i
            if (!found_source && tail->length >= len_lcs - t_i)
            {
                if (source.idx == GVA_NULL)
                {
                    va_array_append(va_std_allocator, graph.nodes, ((Node) {source.row, source.col, source.length, GVA_NULL, GVA_NULL}));
                    source.idx = va_array_length(graph.nodes) - 1;
                    graph.source = source.idx;
                } // if
                printf("SOURCE    %u %u %u: ", source.row, source.col, source.length);
                // is this edges2 check always true?
                if (edges2(source, *tail, true, is_sink, len_obs, observed, &(uint32_t) {0}, &(uint32_t) {0}))
                {
                    va_array_append(va_std_allocator, graph.edges, ((Edge2) {tail->idx, graph.nodes[source.idx].edges}));
                    graph.nodes[source.idx].edges = va_array_length(graph.edges) - 1;
                } // if
            } // if
        } // for t_j
        // dealloc lcs pos level
        lcs_nodes[len_lcs - t_i - 1] = va_array_destroy(va_std_allocator, lcs_nodes[len_lcs - t_i - 1]);
    } // for t_i

    // deallocation
    lcs_nodes = va_std_allocator.alloc(va_std_allocator.context, lcs_nodes, len_lcs, 0);

    printf("#nodes: %zu\n#edges: %zu\n", va_array_length(graph.nodes), va_array_length(graph.edges));
    printf("source: %u\n", graph.source);
    for (size_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        printf("%zu: (%u, %u, %u):\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length);
        if (graph.nodes[i].lambda != GVA_NULL)
        {
            printf("    (%u, %u, %u): lambda\n", graph.nodes[graph.nodes[i].lambda].row, graph.nodes[graph.nodes[i].lambda].col, graph.nodes[graph.nodes[i].lambda].length);
        }
        for (size_t j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
        {
            printf("    (%u, %u, %u): ", graph.nodes[graph.edges[j].tail].row, graph.nodes[graph.edges[j].tail].col, graph.nodes[graph.edges[j].tail].length);
            edges2(((VA_LCS_Node) {graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, 0, 0}),
                   ((VA_LCS_Node) {graph.nodes[graph.edges[j].tail].row, graph.nodes[graph.edges[j].tail].col, graph.nodes[graph.edges[j].tail].length, 0, 0}),
                   graph.nodes[i].row == shift && graph.nodes[i].col == shift, graph.nodes[graph.edges[j].tail].edges == GVA_NULL, len_obs, observed, &(uint32_t) {0}, &(uint32_t) {0});
        } // for
    } // for


    bfs_traversal2(graph, len_obs, observed);

    return graph;
} // build


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
        check(len_ref, reference, len_obs, observed, true);

        build(len_ref, reference, len_obs, observed, 0);

/*
        edges2((VA_LCS_Node) {.row = 1, .col = 5, .length = 5}, (VA_LCS_Node) {.row = 6, .col = 6, .length = 4}, false, true, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 7, .col = 8, .length = 1}, (VA_LCS_Node) {.row = 6, .col = 6, .length = 4}, false, true, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 4, .col = 7, .length = 1}, (VA_LCS_Node) {.row = 6, .col = 6, .length = 4}, false, true, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 1, .length = 1}, (VA_LCS_Node) {.row = 6, .col = 6, .length = 4}, false, true, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 2, .length = 1}, (VA_LCS_Node) {.row = 6, .col = 6, .length = 4}, false, true, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 3, .length = 1}, (VA_LCS_Node) {.row = 6, .col = 6, .length = 4}, false, true, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 4, .length = 1}, (VA_LCS_Node) {.row = 6, .col = 6, .length = 4}, false, true, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 5, .col = 0, .length = 1}, (VA_LCS_Node) {.row = 6, .col = 6, .length = 4}, false, true, len_obs, observed);

        edges2((VA_LCS_Node) {.row = 1, .col = 5, .length = 5}, (VA_LCS_Node) {.row = 7, .col = 8, .length = 1}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 5, .length = 5}, (VA_LCS_Node) {.row = 4, .col = 7, .length = 1}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 1, .length = 1}, (VA_LCS_Node) {.row = 1, .col = 5, .length = 5}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 2, .length = 1}, (VA_LCS_Node) {.row = 1, .col = 5, .length = 5}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 3, .length = 1}, (VA_LCS_Node) {.row = 1, .col = 5, .length = 5}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 4, .length = 1}, (VA_LCS_Node) {.row = 1, .col = 5, .length = 5}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 5, .col = 0, .length = 1}, (VA_LCS_Node) {.row = 1, .col = 5, .length = 5}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 0, .col = 0, .length = 0}, (VA_LCS_Node) {.row = 1, .col = 5, .length = 5}, true, false, len_obs, observed);

        edges2((VA_LCS_Node) {.row = 4, .col = 7, .length = 1}, (VA_LCS_Node) {.row = 7, .col = 8, .length = 1}, false, false, len_obs, observed);

        edges2((VA_LCS_Node) {.row = 1, .col = 1, .length = 1}, (VA_LCS_Node) {.row = 5, .col = 0, .length = 1}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 2, .length = 1}, (VA_LCS_Node) {.row = 5, .col = 0, .length = 1}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 3, .length = 1}, (VA_LCS_Node) {.row = 5, .col = 0, .length = 1}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 4, .length = 1}, (VA_LCS_Node) {.row = 5, .col = 0, .length = 1}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 0, .col = 0, .length = 0}, (VA_LCS_Node) {.row = 5, .col = 0, .length = 1}, true, false, len_obs, observed);

        edges2((VA_LCS_Node) {.row = 1, .col = 1, .length = 1}, (VA_LCS_Node) {.row = 1, .col = 4, .length = 1}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 2, .length = 1}, (VA_LCS_Node) {.row = 1, .col = 4, .length = 1}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 1, .col = 3, .length = 1}, (VA_LCS_Node) {.row = 1, .col = 4, .length = 1}, false, false, len_obs, observed);
        edges2((VA_LCS_Node) {.row = 0, .col = 0, .length = 0}, (VA_LCS_Node) {.row = 1, .col = 4, .length = 1}, true, false, len_obs, observed);
*/

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
