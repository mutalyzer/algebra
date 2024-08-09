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


static void
to_dot(Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    printf("digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\nsi[shape=point,width=.1]\nsi->s%u\n", graph.source);
    for (uint32_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        printf("s%u[label=\"(%u, %u, %u)\"%s]\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, graph.nodes[i].edges == (uint32_t) -1 ? ",peripheries=2" : "");
        if (graph.nodes[i].lambda != (uint32_t) -1)
        {
            printf("s%u->s%u[label=\"&#955;\",style=\"dashed\"]\n", i, graph.nodes[i].lambda);
        } // if
        for (uint32_t j = graph.nodes[i].edges; j != (uint32_t) -1; j = graph.edges[j].next)
        {
            printf("s%u->s%u[label=\"%u:%u/%.*s\"]\n", i, graph.edges[j].tail, graph.edges[j].variant.start, graph.edges[j].variant.end, (int) graph.edges[j].variant.obs_end - graph.edges[j].variant.obs_start, observed + graph.edges[j].variant.obs_start);
        } // for
    } // for
    printf("}\n");
} // to_dot


typedef struct
{
    uint32_t post;
    uint32_t rank;
    uint32_t start;
    uint32_t end;
    uint32_t edge;
    uint32_t next;
} Post_Dom_Table;


static uint32_t
intersection(uint32_t lhs, uint32_t rhs, uint32_t const length, Post_Dom_Table const visited[static length])
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


#define MAX(lhs, rhs) (lhs) > (rhs) ? (lhs) : (rhs)
#define MIN(lhs, rhs) (lhs) < (rhs) ? (lhs) : (rhs)


static void
post_dominators(Graph const graph, uint32_t const head, uint32_t* const rank, uint32_t const start, uint32_t const length, Post_Dom_Table visited[static length])
{
    visited[head].post = -1;
    visited[head].start = start;
    visited[head].end = -1;
    if (graph.nodes[head].lambda != (uint32_t) -1 && visited[graph.nodes[head].lambda].rank == (uint32_t) -1)
    {
        post_dominators(graph, graph.nodes[head].lambda, rank, start, length, visited);
        visited[head].post = graph.nodes[head].lambda;
    } // if

    for (uint32_t i = graph.nodes[head].edges; i != (uint32_t) -1; i = graph.edges[i].next)
    {
        if (visited[graph.edges[i].tail].rank == (uint32_t) -1)
        {
            post_dominators(graph, graph.edges[i].tail, rank, graph.edges[i].variant.end, length, visited);
        } // if
        else
        {
            visited[graph.edges[i].tail].start = MAX(visited[graph.edges[i].tail].start, graph.edges[i].variant.end);
        } // else

        if (visited[head].post == (uint32_t) -1)
        {
            visited[head].post = graph.edges[i].tail;
            visited[head].end = graph.edges[i].variant.start;
        } // if
        else
        {
            visited[head].post = intersection(visited[head].post, graph.edges[i].tail, length, visited);
            visited[head].end = MIN(visited[head].end, graph.edges[i].variant.start);
        } // else
    } // for

    visited[head].rank = *rank;
    *rank += 1;
    printf("visit %u (%u): %d\n", head, visited[head].rank, visited[head].post);
} // post_dominators


static void
local_supremal(VA_Allocator const allocator, Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    uint32_t const length = va_array_length(graph.nodes);
    Post_Dom_Table* visited = allocator.alloc(allocator.context, NULL, 0, sizeof(*visited) * length);
    if (visited == NULL)
    {
        return;
    } // if

    for (uint32_t i = 0; i < length; ++i)
    {
        visited[i].rank = -1;
    } // for

    post_dominators(graph, graph.source, &(uint32_t) {0}, 0, length, visited);

    printf(" # \tpost\trank\tstart\tend\n");
    for (uint32_t i = 0; i < length; ++i)
    {
        printf("%2u:\t%4d\t%4u\t%5u\t%3d\n", i, visited[i].post, visited[i].rank, visited[i].start, visited[i].end);
    } // for

    uint32_t head = graph.source;
    uint32_t shift = 0;
    for (uint32_t tail = visited[head].post; tail != (uint32_t) -1; head = tail, tail = visited[head].post)
    {
        uint32_t const start = visited[head].end;
        uint32_t const end = visited[tail].start;
        VA_Variant const variant = {start, end, graph.nodes[head].col + start - graph.nodes[head].row - shift, graph.nodes[tail].col + end - graph.nodes[tail].row - shift};
        printf("%u:%u/%.*s\n", variant.start, variant.end, (int) variant.obs_end - variant.obs_start, observed + variant.obs_start);
    } // for

    visited = allocator.alloc(allocator.context, visited, sizeof(*visited) * length, 0);
} // local_supremal


static void
local_supremal_it(VA_Allocator const allocator, Graph const graph, size_t const len_obs, char const observed[static len_obs])
{
    uint32_t const length = va_array_length(graph.nodes);
    Post_Dom_Table* visited = allocator.alloc(allocator.context, NULL, 0, sizeof(*visited) * length);
    if (visited == NULL)
    {
        return;
    } // if

    for (uint32_t i = 0; i < length; ++i)
    {
        visited[i].post = -1;
        visited[i].next = -1;
        visited[i].start = 0;
        visited[i].end = -1;
    } // for

    uint32_t rank = 0;
    uint32_t head = graph.source;
    visited[head].edge = graph.nodes[head].edges;

    while (head != (uint32_t) -1)
    {
        //printf("enter %u\n", head);
        uint32_t edge = visited[head].edge;
        if (edge == (uint32_t) -1)
        {
            //printf("visit %u\n", head);
            visited[head].rank = rank;
            rank += 1;
            head = visited[head].next;
            continue;
        } // if

        uint32_t const lambda = graph.nodes[head].lambda;
        if (lambda != (uint32_t) -1 && visited[lambda].next == (uint32_t) -1)
        {
            visited[head].post = intersection(visited[head].post, lambda, length, visited);
            visited[lambda].next = head;
            visited[lambda].edge = graph.nodes[lambda].edges;
            head = lambda;
            continue;
        } // if

        //printf("intersection %d, %u\n", visited[head].post, graph.edges[edge].tail);
        visited[head].post = intersection(visited[head].post, graph.edges[edge].tail, length, visited);

        if (visited[graph.edges[edge].tail].next != (uint32_t) -1)
        {
            visited[head].edge = graph.edges[edge].next;
        } // if
        while (edge != (uint32_t) -1)
        {
            uint32_t const tail = graph.edges[edge].tail;
            visited[head].end = MIN(visited[head].end, graph.edges[edge].variant.start);
            visited[tail].start = MAX(visited[tail].start, graph.edges[edge].variant.end);
            if (visited[tail].next == (uint32_t) -1)
            {
                visited[tail].next = head;
                visited[tail].edge = graph.nodes[tail].edges;
                head = tail;
                break;
            } // if
            visited[head].post = intersection(visited[head].post, tail, length, visited);
            edge = graph.edges[edge].next;
        } // while
    } // while

    printf(" # \tpost\trank\tstart\tend\n");
    for (uint32_t i = 0; i < length; ++i)
    {
        printf("%2u:\t%4d\t%4u\t%5u\t%3d\n", i, visited[i].post, visited[i].rank, visited[i].start, visited[i].end);
    } // for

    head = graph.source;
    uint32_t shift = 0;
    for (uint32_t tail = visited[head].post; tail != (uint32_t) -1; head = tail, tail = visited[head].post)
    {
        uint32_t const start = visited[head].end;
        uint32_t const end = visited[tail].start;
        VA_Variant const variant = {start, end, graph.nodes[head].col + start - graph.nodes[head].row - shift, graph.nodes[tail].col + end - graph.nodes[tail].row - shift};
        printf("%u:%u/%.*s\n", variant.start, variant.end, (int) variant.obs_end - variant.obs_start, observed + variant.obs_start);
    } // for

    visited = allocator.alloc(allocator.context, visited, sizeof(*visited) * length, 0);
} // local_supremal_it


// AAAATA GAAAAGAAA
// AAT ACT
// AGAATTGCTTGAA AGGGTTAAA


typedef struct
{
    uint32_t head;
    uint32_t tail;
    uint32_t length;
    struct
    {
        uint32_t depth;
        uint32_t next;
    } table[];
} Traversal_Generator;


static Traversal_Generator*
traversal_init(VA_Allocator const allocator, Graph const graph)
{
    uint32_t const length = va_array_length(graph.nodes);
    Traversal_Generator* const generator = allocator.alloc(allocator.context, NULL, 0, sizeof(*generator) + sizeof(generator->table[0]) * length);
    if (generator == NULL)
    {
        return NULL;
    } // if

    for (uint32_t i = 0; i < length; ++i)
    {
        generator->table[i].depth = 0;
        generator->table[i].next = -1;
    } // for
    generator->head = graph.source;
    generator->tail = graph.source;
    generator->length = length;
    return generator;
} // traversal_init


static Traversal_Generator*
traversal_destroy(VA_Allocator const allocator, Traversal_Generator* self)
{
    if (self == NULL)
    {
        return NULL;
    } // if
    self = allocator.alloc(allocator.context, self, sizeof(*self) + sizeof(self->table[0]) * self->length, 0);
    return self;
} // traversal_destroy


static uint32_t
traversal_next(Traversal_Generator* const self, Graph const graph)
{
    if (self == NULL)
    {
        return -1;
    } // if

    uint32_t const head = self->head;
    if (head == (uint32_t) -1)
    {
        return head;
    } // if

    for (uint32_t i = graph.nodes[head].edges; i != (uint32_t) -1; i = graph.edges[i].next)
    {
        uint32_t const tail = graph.edges[i].tail;

        uint32_t iets = head;
        for (uint32_t j = graph.nodes[tail].lambda; j != (uint32_t) -1; j = graph.nodes[j].lambda)
        {
            if (self->table[j].depth == 0)
            {
                //printf("push (lambda) %u\n", j);
                self->table[j].depth = self->table[head].depth + 1;
                self->table[j].next = self->table[iets].next;
                self->table[iets].next = j;
                iets = j;
            } // if
        } // for

        if (self->table[tail].depth == 0)
        {
            //printf("push %u\n", tail);
            self->table[tail].depth = self->table[head].depth + 1;
            self->table[tail].next = self->table[head].next;
            self->table[head].next = tail;
            self->head = tail;
        } // if
    } // for

    /*
    printf(" #\tdepth\tnext\n");
    for (uint32_t i = 0; i < self->length; ++i)
    {
        printf("%2u:\t%5u\t%4d\n", i, self->table[i].depth, self->table[i].next);
    } // for
    printf("head: %d\ntail: %d\n", self->head, self->tail);
    */
    return head;
} // traversal_next


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
    printf("digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\nsi[shape=point,width=.1]\n");
    printf("si->s%u\n", graph.source);

    while (head != (uint32_t) -1)
    {
        //printf("pop %u\n", head);
        for (uint32_t i = head; i != (uint32_t) -1; i = graph.nodes[i].lambda)
        {
            printf("s%u[label=\"(%u, %u, %u)\"%s]\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, graph.nodes[i].edges == (uint32_t) -1 ? ",peripheries=2" : "");
            if (i != head)
            {
                //printf("lambda %u\n", i);
            } // if
            for (uint32_t j = graph.nodes[i].edges; j != (uint32_t) -1; j = graph.edges[j].next)
            {
                count += 1;
                printf("s%u->s%u[label=\"%u:%u/%.*s\"]\n", head, graph.edges[j].tail, graph.edges[j].variant.start, graph.edges[j].variant.end, (int) graph.edges[j].variant.obs_end - graph.edges[j].variant.obs_start, observed + graph.edges[j].variant.obs_start);

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

    printf("}\n");

    table = allocator.alloc(allocator.context, table, sizeof(*table) * va_array_length(graph.nodes), 0);

    return count;
} // bfs_traversal


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

    /*
    printf("%zu\n", to_dot(graph, len_obs, observed));
    */

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

    //printf("%zu\n", bfs_traversal(va_std_allocator, graph, len_obs, observed));

    to_dot(graph, len_obs, observed);

/*
    Traversal_Generator* gen = traversal_init(va_std_allocator, graph);
    uint32_t node = traversal_next(gen, graph);
    while (node != (uint32_t) -1)
    {
        printf("yield: %u\n", node);
        node = traversal_next(gen, graph);
    } // if
    gen = traversal_destroy(va_std_allocator, gen);
*/

    local_supremal(va_std_allocator, graph, len_obs, observed);
    local_supremal_it(va_std_allocator, graph, len_obs, observed);

    destroy(va_std_allocator, &graph);

    return EXIT_SUCCESS;
} // main
