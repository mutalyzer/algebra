#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint32_t
#include <stdio.h>      // stderr, fprintf, printf

#include "../include/alloc.h"           // VA_Allocator
#include "../include/array.h"           // va_array_*
#include "../include/variant.h"         // VA_Variant
#include "../include/graph.h"           // Edge, Edge2, Graph, Node

#define MAX(lhs, rhs) (((lhs) > (rhs)) ? (lhs) : (rhs))
#define MIN(lhs, rhs) (((lhs) < (rhs)) ? (lhs) : (rhs))

#define print_variant(variant, observed) variant.start, variant.end, (int) variant.obs_end - variant.obs_start, observed + variant.obs_start
#define VAR_FMT "%u:%u/%.*s"


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