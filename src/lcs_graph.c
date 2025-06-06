#include <inttypes.h>   // intmax_t
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t


#include <stdio.h>      // DEBUG


#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, GVA_Variant, gva_lcs_graph_*
#include "align.h"  // LCS_Alignment, lcs_align
#include "array.h"  // ARRAY_*, array_length
#include "base.h"   // gva_uint, MAX, MIN


gva_uint
gva_edges(GVA_Node const head, GVA_Node const tail,
    bool const is_source, bool const is_sink,
    GVA_Variant* const variant)
{
    intmax_t const row = (intmax_t) head.row - is_source;
    intmax_t const col = (intmax_t) head.col - is_source;
    gva_uint const head_length = head.length + is_source;
    gva_uint const tail_length = tail.length + is_sink;

    intmax_t const offset = MIN((intmax_t) tail.row - row, (intmax_t) tail.col - col) - 1;

    gva_uint const head_offset = offset > 0 ? MIN(head_length, offset + 1) : 1;
    gva_uint const tail_offset = offset < 0 ? MIN(tail_length, -offset) : 0;

    *variant = (GVA_Variant) {
        row + head_offset, tail.row + tail_offset,
        col + head_offset, tail.col + tail_offset
    };
    return MIN(head_length - head_offset, tail_length - tail_offset - 1) + 1;
} // gva_edges


GVA_LCS_Graph
gva_lcs_graph_init(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs],
    size_t const shift)
{
    LCS_Alignment lcs = lcs_align(allocator, len_ref, reference, len_obs, observed);

    GVA_LCS_Graph graph = {.nodes = NULL, .edges = NULL};

    if (lcs.nodes == NULL)
    {
        gva_uint const sink = ARRAY_APPEND(allocator, graph.nodes, ((GVA_Node) {
            len_ref, len_obs, 0, GVA_NULL, GVA_NULL
        })) - 1;
        if (len_ref == 0 && len_obs == 0)
        {
            graph.source = sink;
            return graph;
        } // if

        graph.source = ARRAY_APPEND(allocator, graph.nodes, ((GVA_Node) {
            shift, shift, 0,
            ARRAY_APPEND(allocator, graph.edges, ((GVA_Edge) {sink, GVA_NULL})) - 1,
            GVA_NULL
        })) - 1;
        return graph;
    } // if

    struct
    {
        gva_uint count;
        gva_uint idx;
    }* table = allocator.allocate(allocator.context, NULL, 0, sizeof(*table) * lcs.length);

    for (gva_uint i = 0; i < lcs.length; ++i)
    {
        table[i].count = 0;
        table[i].idx = GVA_NULL;
        fprintf(stderr, "%u: ", i);
        for (gva_uint j = lcs.index[i].head; j != GVA_NULL; j = lcs.nodes[j].next)
        {
            fprintf(stderr, "(%u, %u, %u) ", lcs.nodes[j].row, lcs.nodes[j].col, lcs.nodes[j].length);
        } // for
        fprintf(stderr, "\n");
    } // for

    gva_uint tail_idx = lcs.index[lcs.length - 1].tail;
    LCS_Node sink = lcs.nodes[tail_idx];
    if (sink.row + sink.length == len_ref + shift && sink.col + sink.length == len_obs + shift)
    {
        lcs.nodes[tail_idx].idx = ARRAY_APPEND(allocator, graph.nodes, ((GVA_Node) {
            sink.row, sink.col, sink.length, GVA_NULL, GVA_NULL
        })) - 1;
        lcs.nodes[tail_idx].moved = true;
        sink = lcs.nodes[tail_idx];
    } // if
    else
    {
        sink = (LCS_Node) {.row = len_ref + shift, .col = len_obs + shift, .length = 0};
        sink.idx = ARRAY_APPEND(allocator, graph.nodes, ((GVA_Node) {
            sink.row, sink.col, sink.length, GVA_NULL, GVA_NULL
        })) - 1;
        tail_idx = GVA_NULL;
    } // else
    for (gva_uint i = lcs.index[lcs.length - 1].head; i != tail_idx; i = lcs.nodes[i].next)
    {
        lcs.nodes[i].idx = ARRAY_APPEND(allocator, graph.nodes, ((GVA_Node) {
            lcs.nodes[i].row, lcs.nodes[i].col, lcs.nodes[i].length,
            ARRAY_APPEND(allocator, graph.edges, ((GVA_Edge) {sink.idx, GVA_NULL})) - 1,
            GVA_NULL
        })) - 1;
    } // for

    for (gva_uint i = lcs.length - 1; i >= 1; --i)
    {
        gva_uint next = GVA_NULL;
        for (gva_uint j = lcs.index[i].head; j != GVA_NULL; j = next)
        {
            LCS_Node* const restrict tail = &lcs.nodes[j];
            next = tail->next;
            if (lcs.nodes[j].idx == GVA_NULL)
            {
                continue;
            } // if

            table[i].count += 1;
            table[i].idx = j;

            gva_uint here = GVA_NULL;
            for (gva_uint k = lcs.index[i - 1].head; k != GVA_NULL; k = lcs.nodes[k].next)
            {
                LCS_Node* const restrict head = &lcs.nodes[k];

                if (k >= j ||
                    head->row + head->length >= tail->row + tail->length ||
                    head->col + head->length >= tail->col + tail->length)
                {
                    continue;
                } // if

                here = k;
                if (head->incoming == i)
                {
                    gva_uint const split_idx = head->idx;
                    head->idx = ARRAY_APPEND(allocator, graph.nodes, ((GVA_Node) {
                        head->row, head->col, head->length,
                        ARRAY_APPEND(allocator, graph.edges, ((GVA_Edge) {tail->idx, GVA_NULL})) - 1,
                        split_idx
                    })) - 1;
                    head->moved = false;
                    head->incoming = 0;

                    graph.nodes[split_idx].row += head->length;
                    graph.nodes[split_idx].col += head->length;
                    graph.nodes[split_idx].length -= head->length;
                } // if
                else if (head->idx == GVA_NULL)
                {
                    head->idx = ARRAY_APPEND(allocator, graph.nodes, ((GVA_Node) {
                        head->row, head->col, head->length,
                        ARRAY_APPEND(allocator, graph.edges, ((GVA_Edge) {tail->idx, GVA_NULL})) - 1,
                        GVA_NULL
                    })) - 1;
                } // if
                else if (!head->moved || !tail->moved)
                {
                    graph.nodes[head->idx].edges = ARRAY_APPEND(allocator, graph.edges, ((GVA_Edge) {
                        tail->idx, graph.nodes[head->idx].edges
                    })) - 1;
                } // if
            } // for

            if (tail->length > 1)
            {
                tail->length -= 1;
                tail->moved = true;
                if (here != GVA_NULL)
                {
                    tail->incoming = i;
                    tail->next = lcs.nodes[here].next;
                    lcs.nodes[here].next = j;
                } // if
                else
                {
                    tail->next = lcs.index[i - 1].head;
                    lcs.index[i - 1].head = j;
                } // else
            } // if
        } // for
    } // for

    gva_uint head_idx = lcs.index[0].head;
    LCS_Node source = lcs.nodes[head_idx];
    if (source.row == shift && source.col == shift)
    {
        table[0].count += 1;
        table[0].idx = head_idx;

        head_idx = source.next;
    } // if
    else
    {
        source = (LCS_Node) {.row = shift, .col = shift, .length = 0};
        source.idx = ARRAY_APPEND(allocator, graph.nodes, ((GVA_Node) {
            source.row, source.col, source.length, GVA_NULL, GVA_NULL
        })) - 1;
    } // else
    for (gva_uint i = head_idx; i != GVA_NULL; i = lcs.nodes[i].next)
    {
        if (lcs.nodes[i].idx == GVA_NULL)
        {
            continue;
        } // if

        table[0].count += 1;
        table[0].idx = i;

        if (source.length == 0 || !lcs.nodes[i].moved)
        {
            graph.nodes[source.idx].edges = ARRAY_APPEND(allocator, graph.edges, ((GVA_Edge) {
                lcs.nodes[i].idx, graph.nodes[source.idx].edges
            })) - 1;
        } // if
    } // for

    gva_uint len = 0;
    gva_uint prev = GVA_NULL;
    for (gva_uint i = 0; i < lcs.length; ++i)
    {
        //gva_uint const idx = lcs.nodes[table[i].idx].idx;
        //fprintf(stderr, "%u: %u  (%u, %u, %u)\n", i, table[i].count, graph.nodes[idx].row, graph.nodes[idx].col, graph.nodes[idx].length);
        if ((prev != table[i].idx || table[i].count > 1) && len > 0)
        {
            gva_uint const idx = lcs.nodes[table[i - 1].idx].idx;
            fprintf(stderr, "uniq_match in (%u, %u, %u) @ %u of length %u\n", graph.nodes[idx].row, graph.nodes[idx].col, graph.nodes[idx].length, i - len, len);
            len = 0;
        } // if
        if (table[i].count == 1)
        {
            len += 1;
        } // if
        prev = table[i].idx;
    } // for
    if (len > 0)
    {
        gva_uint const idx = lcs.nodes[table[lcs.length - 1].idx].idx;
        fprintf(stderr, "uniq_match in (%u, %u, %u) @ %zu of length %u\n", graph.nodes[idx].row, graph.nodes[idx].col, graph.nodes[idx].length, lcs.length - len - 1, len);
    } // if

    gva_uint min_source = GVA_NULL;
    for (gva_uint i = graph.nodes[source.idx].edges; i != GVA_NULL; i = graph.edges[i].next)
    {
        GVA_Variant variant;
        gva_edges(graph.nodes[source.idx], graph.nodes[graph.edges[i].tail], true, false, &variant);
        min_source = MIN(min_source, variant.start);
    } // for
    if (min_source != GVA_NULL)
    {
        graph.nodes[source.idx].row += min_source;
        graph.nodes[source.idx].col += min_source;
        graph.nodes[source.idx].length -= min_source;
    } // if
    else
    {
        graph.nodes[source.idx].length = 0;
    } // else

    graph.source = source.idx;

    table = allocator.allocate(allocator.context, table, sizeof(*table) * lcs.length, 0);
    lcs.index = allocator.allocate(allocator.context, lcs.index, lcs.length * sizeof(*lcs.index), 0);
    lcs.nodes = ARRAY_DESTROY(allocator, lcs.nodes);
    return graph;
} // gva_lcs_graph_init


inline void
gva_lcs_graph_destroy(GVA_Allocator const allocator, GVA_LCS_Graph self)
{
    self.nodes = ARRAY_DESTROY(allocator, self.nodes);
    self.edges = ARRAY_DESTROY(allocator, self.edges);
} // gva_lcs_graph_destroy
