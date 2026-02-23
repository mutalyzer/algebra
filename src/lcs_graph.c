#include <inttypes.h>   // intmax_t
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint8_t
#include <string.h>     // memcpy

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_lcs_graph_*
#include "../include/string.h"      // GVA_String, gva_string_*
#include "../include/variant.h"     // GVA_Variant, gva_variant_length
#include "align.h"      // LCS_Alignment, lcs_align
#include "array.h"      // ARRAY_*, array_length
#include "bitset.h"     // bitset_add
#include "common.h"     // GVA_NULL, MAX, MIN, gva_uint


GVA_LCS_Graph
gva_lcs_graph_init(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs],
    size_t const offset)
{
    LCS_Alignment lcs = lcs_align(allocator, len_ref, reference, len_obs, observed, offset);
    size_t const distance = len_ref + len_obs - 2 * lcs.length;

    GVA_LCS_Graph graph =
    {
        .supremal = {0, 0, {0, observed}},
        .observed = {len_obs, observed},
        .source = GVA_NULL,
        .distance = distance,
    };

    if (lcs.nodes == NULL || distance == 0)
    {
        gva_uint const sink = ARRAY_APPEND(allocator, graph.nodes,
            ((GVA_Node)
            {
                .row = len_ref + offset,
                .col = len_obs,
                .edges = GVA_NULL,
                .lambda = GVA_NULL,
            }));

        if (distance == 0)
        {
            graph.source = sink;
            graph.nodes[sink].row = 0;
            graph.nodes[sink].col = 0;
            ARRAY_APPEND(allocator, graph.dom_nodes, ((GVA_Node) {.link = graph.source}));

            lcs.index = allocator.allocate(allocator.context, lcs.index, lcs.length * sizeof(*lcs.index), 0);
            lcs.nodes = ARRAY_DESTROY(allocator, lcs.nodes);
            return graph;
        } // if

        graph.source = ARRAY_APPEND(allocator, graph.nodes,
            ((GVA_Node)
            {
                .row = offset,
                .edges = ARRAY_APPEND(allocator, graph.edges, ((GVA_Edge) {sink, GVA_NULL})),
                .lambda = GVA_NULL,
            }));

        ARRAY_APPEND(allocator, graph.dom_nodes,
             ((GVA_Node)
             {
                 .row = offset,
                 .link = graph.source,
             }));
        ARRAY_APPEND(allocator, graph.dom_nodes,
             ((GVA_Node)
             {
                 .row = len_ref + offset,
                 .col = len_obs,
                 .distance = distance,
                 .link = sink,
             }));
        graph.supremal = (GVA_Variant) {offset, len_ref + offset, {len_obs, observed}};
        return graph;
    } // if

    gva_uint tail_idx = lcs.index[lcs.length - 1].tail;
    LCS_Node sink = lcs.nodes[tail_idx];
    if (sink.row + sink.length == len_ref + offset && sink.col + sink.length == len_obs)
    {
        lcs.nodes[tail_idx].idx = ARRAY_APPEND(allocator, graph.nodes,
            ((GVA_Node)
            {
                .row = sink.row,
                .col = sink.col,
                .length = sink.length,
                .edges = GVA_NULL,
                .lambda = GVA_NULL,
            }));
        lcs.nodes[tail_idx].moved = true;
        sink = lcs.nodes[tail_idx];
    } // if
    else
    {
        sink = (LCS_Node) {.row = len_ref + offset, .col = len_obs, .length = 0};
        sink.idx = ARRAY_APPEND(allocator, graph.nodes,
            ((GVA_Node)
            {
                .row = sink.row,
                .col = sink.col,
                .length = sink.length,
                .edges = GVA_NULL,
                .lambda = GVA_NULL,
            }));
        tail_idx = GVA_NULL;
    } // else
    for (gva_uint i = lcs.index[lcs.length - 1].head; i != tail_idx; i = lcs.nodes[i].next)
    {
        lcs.nodes[i].idx = ARRAY_APPEND(allocator, graph.nodes,
            ((GVA_Node)
            {
                .row = lcs.nodes[i].row,
                .col = lcs.nodes[i].col,
                .length = lcs.nodes[i].length,
                .edges = ARRAY_APPEND(allocator, graph.edges, ((GVA_Edge) {sink.idx, GVA_NULL})),
                .lambda = GVA_NULL,
            }));
    } // for

    for (gva_uint i = lcs.length - 1; i >= 1; --i)
    {
        gva_uint next = GVA_NULL;
        for (gva_uint j = lcs.index[i].head; j != GVA_NULL; j = next)
        {
            LCS_Node* const restrict tail = &lcs.nodes[j];
            next = tail->next;
            if (tail->idx == GVA_NULL)
            {
                continue;
            } // if

            lcs.index[i].tail = j;
            lcs.index[i].count += 1;
            lcs.index[i].offset = graph.nodes[tail->idx].length - tail->length;

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
                    head->idx = ARRAY_APPEND(allocator, graph.nodes,
                        ((GVA_Node)
                        {
                            .row = head->row,
                            .col = head->col,
                            .length = head->length,
                            .edges = ARRAY_APPEND(allocator, graph.edges, ((GVA_Edge) {tail->idx, GVA_NULL})),
                            .lambda = split_idx
                        }));
                    head->moved = false;
                    head->incoming = 0;

                    graph.nodes[split_idx].row += head->length;
                    graph.nodes[split_idx].col += head->length;
                    graph.nodes[split_idx].length -= head->length;
                } // if
                else if (head->idx == GVA_NULL)
                {
                    head->idx = ARRAY_APPEND(allocator, graph.nodes,
                        ((GVA_Node)
                        {
                            .row = head->row,
                            .col = head->col,
                            .length = head->length,
                            .edges = ARRAY_APPEND(allocator, graph.edges, ((GVA_Edge) {tail->idx, GVA_NULL})),
                            .lambda = GVA_NULL
                        }));
                } // if
                else if (!head->moved || !tail->moved)
                {
                    graph.nodes[head->idx].edges = ARRAY_APPEND(allocator, graph.edges,
                        ((GVA_Edge)
                        {
                            .tail = tail->idx,
                            .next = graph.nodes[head->idx].edges,
                        }));
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
    if (source.row == offset && source.col == 0)
    {
        lcs.index[0].tail = head_idx;
        lcs.index[0].count += 1;
        lcs.index[0].offset = graph.nodes[source.idx].length - source.length;

        head_idx = source.next;
    } // if
    else
    {
        source = (LCS_Node) {.row = offset, .col = 0, .length = 0};
        source.idx = ARRAY_APPEND(allocator, graph.nodes,
            ((GVA_Node)
            {
                .row = source.row,
                .col = source.col,
                .length = source.length,
                .edges = GVA_NULL,
                .lambda = GVA_NULL,
            }));
    } // else
    for (gva_uint i = head_idx; i != GVA_NULL; i = lcs.nodes[i].next)
    {
        if (lcs.nodes[i].idx == GVA_NULL)
        {
            continue;
        } // if

        lcs.index[0].tail = i;
        lcs.index[0].count += 1;
        lcs.index[0].offset = graph.nodes[lcs.nodes[i].idx].length - lcs.nodes[i].length;

        if (source.length == 0 || !lcs.nodes[i].moved)
        {
            graph.nodes[source.idx].edges = ARRAY_APPEND(allocator, graph.edges,
                ((GVA_Edge)
                {
                    .tail = lcs.nodes[i].idx,
                    .next = graph.nodes[source.idx].edges,
                }));
        } // if
    } // for

    // start constructing the local supremal
    ARRAY_APPEND(allocator, graph.dom_nodes,
        ((GVA_Node)
        {
            .row = offset,
            .link = source.idx,
        }));
    gva_uint prev = lcs.index[0].tail;
    if (lcs.nodes[lcs.index[0].tail].idx != source.idx)
    {
        prev = GVA_NULL;
    } // if
    for (gva_uint i = 0; i < lcs.length; ++i)
    {
        // when to extend
        if (lcs.index[i].count == 1)
        {
            // start new
            if (lcs.index[i].tail != prev)
            {
                gva_uint const pos = graph.nodes[lcs.nodes[lcs.index[i].tail].idx].length - lcs.index[i].offset - 1;
                ARRAY_APPEND(allocator, graph.dom_nodes,
                    ((GVA_Node)
                    {
                        .row = lcs.nodes[lcs.index[i].tail].row + pos,
                        .col = lcs.nodes[lcs.index[i].tail].col + pos,
                        .distance = lcs.nodes[lcs.index[i].tail].row + pos - offset +
                                    lcs.nodes[lcs.index[i].tail].col + pos - 2 * i,
                        .link = lcs.nodes[lcs.index[i].tail].idx,
                    }));

                prev = lcs.index[i].tail;
            } // if

            graph.dom_nodes[array_length(graph.dom_nodes) - 1].length += 1;
        } // if
    } // for

    // add empty match sink if not already handled because of unique matches
    if (graph.dom_nodes[array_length(graph.dom_nodes) - 1].link != sink.idx)
    {
        ARRAY_APPEND(allocator, graph.dom_nodes,
            ((GVA_Node)
            {
                .row = sink.row + sink.length,
                .col = sink.col + sink.length,
                .distance = distance,
                .link = sink.idx,
            }));
    } // if

    // construct supremal
    if (graph.dom_nodes != NULL)
    {
        graph.nodes[source.idx].row += graph.dom_nodes[0].length;
        graph.nodes[source.idx].col += graph.dom_nodes[0].length;
        graph.nodes[source.idx].length -= graph.dom_nodes[0].length;

        graph.nodes[sink.idx].length -= graph.dom_nodes[array_length(graph.dom_nodes) - 1].length;

        gva_edges(graph.observed.str,
            graph.dom_nodes[0], graph.dom_nodes[array_length(graph.dom_nodes) - 1],
            true, true,
            &graph.supremal);
    } // if

    graph.source = source.idx;

    lcs.index = allocator.allocate(allocator.context, lcs.index, lcs.length * sizeof(*lcs.index), 0);
    lcs.nodes = ARRAY_DESTROY(allocator, lcs.nodes);
    return graph;
} // gva_lcs_graph_init


inline void
gva_lcs_graph_destroy(GVA_Allocator const allocator, GVA_LCS_Graph self)
{
    self.nodes = ARRAY_DESTROY(allocator, self.nodes);
    self.edges = ARRAY_DESTROY(allocator, self.edges);
    self.dom_nodes = ARRAY_DESTROY(allocator, self.dom_nodes);
    // FIXME: ownership of observed
} // gva_lcs_graph_destroy


GVA_LCS_Graph
gva_lcs_graph_from_variants(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const n, GVA_Variant const variants[static restrict n])
{
    GVA_Variant variant = {variants[0].start, variants[n - 1].end, {0, NULL}};
    variant.sequence = gva_string_concat(allocator, variant.sequence, variants[0].sequence);
    for (size_t i = 1; i < n; ++i)
    {
        variant.sequence = gva_string_concat(allocator, variant.sequence, (GVA_String) {variants[i].start - variants[i - 1].end, reference + variants[i - 1].end});
        variant.sequence = gva_string_concat(allocator, variant.sequence, variants[i].sequence);
    } // for

    gva_uint offset = MAX(8, gva_variant_length(variant) / 2);
    size_t old_len = 0;
    char* observed = NULL;
    while (true)
    {
        gva_uint const start = MAX(0, (intmax_t) variant.start - offset);
        gva_uint const end = MIN(len_ref, variant.end + offset);

        // FIXME: len == 0?
        size_t const len = (variant.start - start) + variant.sequence.len + (end - variant.end);
        observed = allocator.allocate(allocator.context, observed, old_len, len);
        if (observed == NULL)
        {
            gva_string_destroy(allocator, variant.sequence);
            return (GVA_LCS_Graph) {NULL, NULL, NULL, {0, 0, {0, NULL}}, {0, NULL}, GVA_NULL, 0};
        } // if

        memcpy(observed, reference + start, variant.start - start);
        memcpy(observed + variant.start - start, variant.sequence.str, variant.sequence.len);
        memcpy(observed + variant.start - start + variant.sequence.len, reference + variant.end, end - variant.end);

        GVA_LCS_Graph graph = gva_lcs_graph_init(allocator, end - start, reference + start, len, observed, start);

        if ((graph.supremal.start > start || graph.supremal.start == 0) &&
            (graph.supremal.end   < end   || graph.supremal.end   == len_ref))
        {
            // FIXME: observed is now owned by graph
            gva_string_destroy(allocator, variant.sequence);
            return graph;
        } // if

        gva_lcs_graph_destroy(allocator, graph);

        old_len = len;
        offset *= 2;  // OVERFLOW
    } // while
} // gva_lcs_graph_from_variants


static uint8_t const NUC_A = 0x1;
static uint8_t const NUC_C = 0x2;
static uint8_t const NUC_G = 0x4;
static uint8_t const NUC_T = 0x8;


inline static uint8_t
nucleotides(size_t const len, char const sequence[static len])
{
    static uint8_t const MASK[256] =
    {
        ['A'] = NUC_A,
        ['C'] = NUC_C,
        ['G'] = NUC_G,
        ['T'] = NUC_T,
    };
    static uint8_t const UNIVERSE = 0xF;

    uint8_t mask = 0x0;
    for (size_t i = 0; mask < UNIVERSE && i < len; ++i)
    {
        mask |= MASK[(size_t) sequence[i]];
    } // for
    return mask;
} // nucleotides


void
gva_lcs_graph_uniq_atomics(GVA_LCS_Graph const self,
    gva_uint const offset,
    gva_uint const start, gva_uint const end,
    size_t dels[static restrict 1],
    size_t as[static restrict 1],
    size_t cs[static restrict 1],
    size_t gs[static restrict 1],
    size_t ts[static restrict 1])
{
    for (size_t i = 0; i < array_length(self.nodes); ++i)
    {
        if (self.nodes[i].row > end)
        {
            continue;
        } // if

        for (gva_uint j = self.nodes[i].edges; j != GVA_NULL; j = self.edges[j].next)
        {
            if (self.nodes[self.edges[j].tail].row + self.nodes[self.edges[j].tail].length < start)
            {
                continue;
            } // if

            GVA_Variant variant;
            gva_uint const count = gva_edges(self.observed.str,
                                             self.nodes[i], self.nodes[self.edges[j].tail],
                                             i == self.source, self.nodes[self.edges[j].tail].edges == GVA_NULL,
                                             &variant);

            if (variant.end > variant.start)
            {
                bitset_add(dels, variant.start - offset, variant.end + count - 1 - offset);
            } // if

            uint8_t const mask = nucleotides(variant.sequence.len, variant.sequence.str);
            if ((mask & NUC_A) == NUC_A)
            {
                bitset_add(as, variant.start - offset, variant.end + count - offset);
            } // if
            if ((mask & NUC_C) == NUC_C)
            {
                bitset_add(cs, variant.start - offset, variant.end + count - offset);
            } // if
            if ((mask & NUC_G) == NUC_G)
            {
                bitset_add(gs, variant.start - offset, variant.end + count - offset);
            } // if
            if ((mask & NUC_T) == NUC_T)
            {
                bitset_add(ts, variant.start - offset, variant.end + count - offset);
            } // if
        } // for
    } // for
} // gva_lcs_graph_uniq_atomics


inline GVA_Variant
gva_lcs_graph_ls_slice(GVA_LCS_Graph const self,
    size_t const start, size_t const end)
{
    GVA_Variant variant;
    gva_edges(self.observed.str,
              self.dom_nodes[start], self.dom_nodes[end],
              start == 0, end == array_length(self.dom_nodes) - 1,
              &variant);
    return variant;
} // gva_lcs_graph_ls_slice


inline size_t
gva_lcs_graph_distance(GVA_LCS_Graph const self)
{
    return self.dom_nodes[array_length(self.dom_nodes) - 1].distance;
} // gva_lcs_graph_distance


inline GVA_Variant
gva_lcs_graph_supremal(GVA_LCS_Graph const self)
{
    GVA_Variant variant = {0, 0, {0, NULL}};
    gva_edges(self.observed.str,
        self.dom_nodes[0], self.dom_nodes[array_length(self.dom_nodes) - 1],
        true, true,
        &variant);
    return variant;
} // gva_lcs_graph_supremal


gva_uint
gva_edges(char const observed[static restrict 1],
    GVA_Node const head, GVA_Node const tail,
    bool const is_source, bool const is_sink,
    GVA_Variant variant[static restrict 1])
{
    intmax_t const row = (intmax_t) head.row - is_source;
    intmax_t const col = (intmax_t) head.col - is_source;
    gva_uint const head_length = head.length + is_source;
    gva_uint const tail_length = tail.length + is_sink;

    intmax_t const offset = MIN((intmax_t) tail.row - row, (intmax_t) tail.col - col) - 1;

    gva_uint const head_offset = offset > 0 ? MIN(head_length, offset + 1) : 1;
    gva_uint const tail_offset = offset < 0 ? MIN(tail_length, -offset) : 0;

    *variant = (GVA_Variant)
    {
        .start = row + head_offset,
        .end = tail.row + tail_offset,
        .sequence = {(tail.col + tail_offset) - (col + head_offset), observed + col + head_offset},
    };
    return MIN(head_length - head_offset, tail_length - tail_offset - 1) + 1;
} // gva_edges
