#include <stddef.h>     // NULL, size_t


#include <stdio.h>              // DEBUG
#include "../include/utils.h"   // DEBUG


#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_edges
#include "array.h"      // array_length
#include "common.h"     // GVA_NULL, MAX, MIN, GVA_String, gva_uint


typedef struct
{
    gva_uint lca;
    gva_uint rank;
    gva_uint depth;
    gva_uint start;
    gva_uint end;
    gva_uint prev;
    gva_uint next;
} LCA_Table;


static void
print_lca_table(size_t const length, LCA_Table const visited[static length], gva_uint const head, gva_uint const tail)
{
    fprintf(stderr, " # \tlca\trank\tdepth\tstart\tend\tprev\tnext\n");
    for (gva_uint i = 0; i < length; ++i)
    {
        fprintf(stderr, "%2u:\t%3d\t%4d\t%5d\t%5d\t%3d\t%4d\t%4d\n", i, visited[i].lca, visited[i].rank, visited[i].depth, visited[i].start, visited[i].end, visited[i].prev, visited[i].next);
    } // for
    fprintf(stderr, "head: %d\ntail: %d\n", head, tail);
} // print_lca_table


static inline gva_uint
lca(gva_uint* const start, gva_uint lhs, gva_uint rhs, size_t const length, LCA_Table const visited[static length])
{
    gva_uint lhs_start = *start;
    gva_uint rhs_start = *start;
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


void
gva_extract(GVA_Allocator const allocator, GVA_LCS_Graph const graph)
{
    gva_uint const length = array_length(graph.nodes);
    LCA_Table* visited = allocator.allocate(allocator.context, NULL, 0, sizeof(*visited) * length);
    if (visited == NULL)
    {
        return;
    } // if

    for (gva_uint i = 0; i < length; ++i)
    {
        visited[i].depth = GVA_NULL;
        visited[i].lca = GVA_NULL;
        visited[i].rank = GVA_NULL;
        visited[i].start = GVA_NULL;
        visited[i].end = GVA_NULL;
        visited[i].prev = GVA_NULL;
        visited[i].next = GVA_NULL;
    } // for

    gva_uint sink = GVA_NULL;
    visited[graph.source] = (LCA_Table) {GVA_NULL, 0, 0, GVA_NULL, GVA_NULL, GVA_NULL, GVA_NULL};
    gva_uint rank = 1;
    // main loop over a queue
    for (gva_uint head = graph.source, tail = graph.source; head != GVA_NULL; head = visited[head].next)
    {
        fprintf(stderr, "pop %u @ %u\n", head, visited[head].depth);

        if (graph.nodes[head].edges == GVA_NULL)
        {
            fprintf(stderr, "    sink\n");
            sink = head;
            continue;
        } // if

        gva_uint const lambda = graph.nodes[head].lambda;
        if (lambda != GVA_NULL)
        {
            if (visited[lambda].depth == GVA_NULL)
            {
                // add lambda in stack order
                fprintf(stderr, "    push lambda %u @ %u\n", lambda, visited[head].depth);

                visited[visited[head].next].prev = lambda;
                visited[lambda] = (LCA_Table) {head, rank, visited[head].depth, GVA_NULL, GVA_NULL, head, visited[head].next};
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
                visited[lambda] = (LCA_Table) {head, visited[lambda].rank, visited[head].depth, GVA_NULL, GVA_NULL, visited[lambda].prev, visited[lambda].next};
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
                visited[lambda] = (LCA_Table) {head, visited[lambda].rank, visited[head].depth, GVA_NULL, GVA_NULL, head, visited[head].next};
                visited[head].next = lambda;
            } // else

            //print_lca_table(length, visited, head, tail);
        } // if

        for (gva_uint i = graph.nodes[head].edges; i != GVA_NULL; i = graph.edges[i].next)
        {
            gva_uint const edge_tail = graph.edges[i].tail;
            if (visited[edge_tail].depth == GVA_NULL)
            {
                // add regular successors in queue order
                fprintf(stderr, "    push %u @ %u\n", edge_tail, visited[head].depth + 1);
                GVA_Variant variant;
                gva_uint const count = gva_edges(graph.observed.str, graph.nodes[head], graph.nodes[edge_tail], head == graph.source, edge_tail == 0, &variant);
                visited[edge_tail] = (LCA_Table) {head, rank, visited[head].depth + 1, variant.start, variant.end + count - 1, tail, GVA_NULL};
                rank += 1;
                visited[tail].next = edge_tail;
                tail = edge_tail;
            } // if
            else if (visited[edge_tail].depth - 1 == visited[head].depth)
            {
                fprintf(stderr, "    update %u @ %u\n", edge_tail, visited[edge_tail].depth);
                fprintf(stderr, "        lca(%u, %u)\n", visited[edge_tail].lca, head);
                GVA_Variant variant;
                gva_uint const count = gva_edges(graph.observed.str, graph.nodes[head], graph.nodes[edge_tail], head == graph.source, edge_tail == 0, &variant);
                visited[edge_tail].start = MIN(visited[edge_tail].start, variant.start);
                visited[edge_tail].lca = lca(&visited[edge_tail].start, visited[edge_tail].lca, head, length, visited);
                visited[edge_tail].end = MAX(visited[edge_tail].end, variant.end + count - 1);
                fprintf(stderr, "        = %u (%u:%u)\n", visited[edge_tail].lca, visited[edge_tail].start, visited[edge_tail].end);
            } // if
            else
            {
                fprintf(stderr, "    skip %u @ %u\n", edge_tail, visited[edge_tail].depth);
            } // else

            //print_lca_table(length, visited, head, tail);
        } // for
    } // for

    print_lca_table(length, visited, GVA_NULL, GVA_NULL);

    // the canonical variant is given in the reverse order (for now): the Python checker needs to reverse this list
    for (gva_uint tail = sink; tail != GVA_NULL; tail = visited[tail].lca)
    {
        gva_uint const head = visited[tail].lca;
        if (head != GVA_NULL && visited[tail].start != GVA_NULL)
        {
            // we need to skip lambda edges: start && end == -1
            gva_uint const start_offset = visited[tail].start - graph.nodes[head].row;
            gva_uint const end_offset = visited[tail].end - graph.nodes[tail].row;
            GVA_Variant const variant = {
                visited[tail].start, visited[tail].end, {
                    graph.observed.str + graph.nodes[head].col + start_offset,
                    (graph.nodes[tail].col + end_offset) - (graph.nodes[head].col + start_offset)
                }
            };

            if (tail != sink)
            {
                fprintf(stderr, ",\n");
            } // if
            fprintf(stderr, "        \"" GVA_VARIANT_FMT "\"", GVA_VARIANT_PRINT(variant));
        } // if
    } // while

    visited = allocator.allocate(allocator.context, visited, sizeof(*visited) * length, 0);
} // gva_extract
