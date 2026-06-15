#include <inttypes.h>   // intmax_t
#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint8_t
#include <string.h>     // memcpy

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, GVA_Dom_Node,
                                    // GVA_Edge, GVA_Node, gva_lcs_graph_*
#include "../include/string.h"      // GVA_String, gva_string_*
#include "../include/variant.h"     // GVA_Variant, gva_variant_length
#include "align.h"      // LCS_Alignment, LCS_Matches, lcs_align*
#include "array.h"      // ARRAY_*, array_length
#include "common.h"     // GVA_NULL, MAX, MIN, gva_uint


#include <stdio.h>      // FIXME: DEBUG


static void
concat(GVA_Allocator const allocator,
    GVA_LCS_Graph lhs[static 1], GVA_LCS_Graph const rhs,
    size_t const offset)
{
    size_t const offset_nodes = array_length(lhs->nodes);
    size_t const offset_edges = array_length(lhs->edges);

    gva_uint sink_idx = GVA_NULL;
    size_t offset_distance = 0;
    if (offset_nodes > 0)
    {
        sink_idx = lhs->dom_nodes[array_length(lhs->dom_nodes) - 1].link;
        offset_distance = lhs->dom_nodes[array_length(lhs->dom_nodes) - 1].distance;

        // merge the sink of lhs with the source of rhs
        lhs->nodes[sink_idx].match.length += (rhs.nodes[rhs.source].match.row + rhs.nodes[rhs.source].match.length) -
            (lhs->nodes[sink_idx].match.row + lhs->nodes[sink_idx].match.length);
        lhs->nodes[sink_idx].edges = rhs.nodes[rhs.source].edges + offset_edges;
        lhs->dom_nodes[array_length(lhs->dom_nodes) - 1].match.length += (rhs.dom_nodes[0].match.row + rhs.dom_nodes[0].match.length) -
            (lhs->dom_nodes[array_length(lhs->dom_nodes) - 1].match.row + lhs->dom_nodes[array_length(lhs->dom_nodes) - 1].match.length);
    } // if
    else
    {
        lhs->source = rhs.source;
    } // else

    for (size_t i = 0; i < array_length(rhs.nodes); ++i)
    {
        if (i != rhs.source || sink_idx == GVA_NULL)
        {
            ARRAY_APPEND(allocator, lhs->nodes,
            ((GVA_Node)
            {
                .match = {rhs.nodes[i].match.row, rhs.nodes[i].match.col + offset, rhs.nodes[i].match.length},
                .edges = rhs.nodes[i].edges == GVA_NULL ? GVA_NULL : rhs.nodes[i].edges + offset_edges,
                .lambda = rhs.nodes[i].lambda == GVA_NULL ? GVA_NULL : rhs.nodes[i].lambda + offset_nodes - (rhs.nodes[i].lambda > rhs.source && sink_idx != GVA_NULL),
            }));
        } // if
    } // for

    for (size_t i = 0; i < array_length(rhs.edges); ++i)
    {
        ARRAY_APPEND(allocator, lhs->edges,
        ((GVA_Edge)
        {
            .tail = rhs.edges[i].tail + offset_nodes - (rhs.edges[i].tail > rhs.source && sink_idx != GVA_NULL),
            .next = rhs.edges[i].next == GVA_NULL ? GVA_NULL : rhs.edges[i].next + offset_edges,
        }));
    } // for

    for (size_t i = (sink_idx != GVA_NULL); i < array_length(rhs.dom_nodes); ++i)
    {
        ARRAY_APPEND(allocator, lhs->dom_nodes,
        ((GVA_Dom_Node)
        {
            .match = {rhs.dom_nodes[i].match.row, rhs.dom_nodes[i].match.col + offset, rhs.dom_nodes[i].match.length},
            .distance = rhs.dom_nodes[i].distance + offset_distance,
            .link = rhs.dom_nodes[i].link + offset_nodes  - (rhs.dom_nodes[i].link > rhs.source && sink_idx != GVA_NULL),
        }));
    } // for
} // concat


/*
int multi;
int recurse;
int base;
*/

// FIXME: recursion!
static void
local_supremal(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs],
    size_t const offset_row, size_t const offset_col,
    GVA_LCS_Graph graph[static restrict 1])
{
    fprintf(stderr, GVA_STRING_FMT " vs " GVA_STRING_FMT "\n", GVA_STRING_PRINT(((GVA_String) {len_ref, reference})), GVA_STRING_PRINT(((GVA_String) {len_obs, observed})));

    if (len_ref == 0 || len_obs == 0)
    {
        //base += 1;
        GVA_LCS_Graph local = gva_lcs_graph_init(allocator, len_ref, reference, len_obs, observed, offset_row);
        fprintf(stderr, "BASE: " GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(gva_lcs_graph_supremal(local)));
        concat(allocator, graph, local, offset_col);
        gva_lcs_graph_destroy(allocator, local, false);
        return;
    } // if

    fprintf(stderr, "F\n");
    LCS_Matches forward = lcs_align_one(allocator, len_ref, reference, len_obs, observed);
    gva_string_reverse(len_ref, (char*) reference);
    gva_string_reverse(len_obs, (char*) observed);

    fprintf(stderr, "B\n");
    LCS_Matches backward = lcs_align_one(allocator, len_ref, reference, len_obs, observed);
    gva_string_reverse(len_ref, (char*) reference);
    gva_string_reverse(len_obs, (char*) observed);

    size_t sum = 0;
    size_t prev_row = -1;
    size_t prev_col = -1;
    for (size_t i = 0; i < forward.max_lcs_pos; ++i)
    {
        size_t const j = forward.max_lcs_pos - i - 1;

        fprintf(stderr, "%zu: (%u, %u) %d    (%zu, %zu) %d\n", i,
            forward.match[i].row, forward.match[i].col, forward.uniq[i],
            len_ref - backward.match[j].row - 1, len_obs - backward.match[j].col - 1, backward.uniq[j]);

        if (forward.match[i].row == len_ref - backward.match[j].row - 1 &&
            forward.match[i].col == len_obs - backward.match[j].col - 1 // &&
            // (forward.uniq[i] == 1 || backward.uniq[j] == 1)
           )
        {
            size_t const distance = forward.match[i].row + forward.match[i].col - 2 * i - sum;
            if (distance > 0)
            {
                //recurse += 1;
                sum += distance;
                fprintf(stderr, "RECURSE %u %u %d %d {\n", forward.match[i].row, forward.match[i].col, forward.uniq[i], backward.uniq[j]);
                local_supremal(allocator,
                    forward.match[i].row - prev_row - 1, reference + prev_row + 1,
                    forward.match[i].col - prev_col - 1, observed + prev_col + 1,
                    offset_row + prev_row + 1, offset_col + prev_col + 1,
                    graph);
                fprintf(stderr, "}\n");
            } // if

            prev_row = forward.match[i].row;
            prev_col = forward.match[i].col;
        } // if
    } // for
    size_t const distance = len_ref + len_obs - 2 * forward.max_lcs_pos - sum;
    if (distance > 0)
    {
        if (prev_row != (size_t) -1)
        {
            //recurse += 1;
            fprintf(stderr, "RECURSE {\n");
            local_supremal(allocator,
                len_ref - prev_row - 1, reference + prev_row + 1,
                len_obs - prev_col - 1, observed + prev_col + 1,
                offset_row + prev_row + 1, offset_col + prev_col + 1,
                graph);
            fprintf(stderr, "}\n");
        } // if
        else
        {
            GVA_LCS_Graph local = gva_lcs_graph_init(allocator, len_ref - prev_row - 1, reference + prev_row + 1, len_obs - prev_col - 1, observed + prev_col + 1, offset_row + prev_row + 1);
            fprintf(stderr, GVA_VARIANT_FMT " %zu\n", GVA_VARIANT_PRINT(gva_lcs_graph_supremal(local)), gva_variant_length(gva_lcs_graph_supremal(local)));
            /*
            if (array_length(local.dom_nodes) == 3)
            {
                multi = 1;
                fprintf(stderr, "MULTI\n");
                for (size_t i = 0; i < array_length(local.dom_nodes) - 1; ++i)
                {
                    fprintf(stderr, "  " GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(gva_lcs_graph_local_supremal(local, i, i + 1)));
                } // for
            } // if
            */
            concat(allocator, graph, local, offset_col + prev_col + 1);
            gva_lcs_graph_destroy(allocator, local, false);
        } // else
    } // if

    backward.match = allocator.allocate(allocator.context, backward.match, MIN(len_ref, len_obs), 0);
    backward.uniq = allocator.allocate(allocator.context, backward.uniq, MIN(len_ref, len_obs), 0);
    forward.match = allocator.allocate(allocator.context, forward.match, MIN(len_ref, len_obs), 0);
    forward.uniq = allocator.allocate(allocator.context, forward.uniq, MIN(len_ref, len_obs), 0);
} // local_supremal


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
        .observed = {len_obs, observed},
        .source = GVA_NULL,
    };

    if (lcs.nodes == NULL || distance == 0)
    {
        gva_uint const sink = ARRAY_APPEND(allocator, graph.nodes,
            ((GVA_Node)
            {
                .match = {len_ref + offset, len_obs},
                .edges = GVA_NULL,
                .lambda = GVA_NULL,
            }));

        if (distance == 0)
        {
            graph.source = sink;
            graph.nodes[sink].match.row = 0;
            graph.nodes[sink].match.col = 0;
            ARRAY_APPEND(allocator, graph.dom_nodes, ((GVA_Dom_Node) {.link = graph.source}));

            lcs.index = allocator.allocate(allocator.context, lcs.index, lcs.length * sizeof(*lcs.index), 0);
            lcs.nodes = ARRAY_DESTROY(allocator, lcs.nodes);
            return graph;
        } // if

        graph.source = ARRAY_APPEND(allocator, graph.nodes,
            ((GVA_Node)
            {
                .match.row = offset,
                .edges = ARRAY_APPEND(allocator, graph.edges, ((GVA_Edge) {sink, GVA_NULL})),
                .lambda = GVA_NULL,
            }));

        ARRAY_APPEND(allocator, graph.dom_nodes,
             ((GVA_Dom_Node)
             {
                 .match.row = offset,
                 .link = graph.source,
             }));
        ARRAY_APPEND(allocator, graph.dom_nodes,
             ((GVA_Dom_Node)
             {
                 .match = {len_ref + offset, len_obs},
                 .distance = distance,
                 .link = sink,
             }));
        return graph;
    } // if

    gva_uint tail_idx = lcs.index[lcs.length - 1].tail;
    LCS_Node sink = lcs.nodes[tail_idx];
    if (sink.match.row + sink.match.length == len_ref + offset &&
        sink.match.col + sink.match.length == len_obs)
    {
        lcs.nodes[tail_idx].idx = ARRAY_APPEND(allocator, graph.nodes,
            ((GVA_Node)
            {
                .match = sink.match,
                .edges = GVA_NULL,
                .lambda = GVA_NULL,
            }));
        lcs.nodes[tail_idx].moved = true;
        sink = lcs.nodes[tail_idx];
    } // if
    else
    {
        sink = (LCS_Node) {.match = {len_ref + offset, len_obs}};
        sink.idx = ARRAY_APPEND(allocator, graph.nodes,
            ((GVA_Node)
            {
                .match = sink.match,
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
                .match = lcs.nodes[i].match,
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
            lcs.index[i].offset = graph.nodes[tail->idx].match.length - tail->match.length;

            gva_uint here = GVA_NULL;
            for (gva_uint k = lcs.index[i - 1].head; k != GVA_NULL; k = lcs.nodes[k].next)
            {
                LCS_Node* const restrict head = &lcs.nodes[k];

                if (k >= j ||
                    head->match.row + head->match.length >= tail->match.row + tail->match.length ||
                    head->match.col + head->match.length >= tail->match.col + tail->match.length)
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
                            .match = head->match,
                            .edges = ARRAY_APPEND(allocator, graph.edges, ((GVA_Edge) {tail->idx, GVA_NULL})),
                            .lambda = split_idx
                        }));
                    head->moved = false;
                    head->incoming = 0;

                    graph.nodes[split_idx].match.row += head->match.length;
                    graph.nodes[split_idx].match.col += head->match.length;
                    graph.nodes[split_idx].match.length -= head->match.length;
                } // if
                else if (head->idx == GVA_NULL)
                {
                    head->idx = ARRAY_APPEND(allocator, graph.nodes,
                        ((GVA_Node)
                        {
                            .match = head->match,
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

            if (tail->match.length > 1)
            {
                tail->match.length -= 1;
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
    if (source.match.row == offset && source.match.col == 0)
    {
        lcs.index[0].tail = head_idx;
        lcs.index[0].count += 1;
        lcs.index[0].offset = graph.nodes[source.idx].match.length - source.match.length;

        head_idx = source.next;
    } // if
    else
    {
        source = (LCS_Node) {.match.row = offset};
        source.idx = ARRAY_APPEND(allocator, graph.nodes,
            ((GVA_Node)
            {
                .match = source.match,
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
        lcs.index[0].offset = graph.nodes[lcs.nodes[i].idx].match.length - lcs.nodes[i].match.length;

        if (source.match.length == 0 || !lcs.nodes[i].moved)
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
        ((GVA_Dom_Node)
        {
            .match.row = offset,
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
                gva_uint const pos = graph.nodes[lcs.nodes[lcs.index[i].tail].idx].match.length - lcs.index[i].offset - 1;
                ARRAY_APPEND(allocator, graph.dom_nodes,
                    ((GVA_Dom_Node)
                    {
                        .match = {lcs.nodes[lcs.index[i].tail].match.row + pos, lcs.nodes[lcs.index[i].tail].match.col + pos},
                        .distance = lcs.nodes[lcs.index[i].tail].match.row + pos - offset +
                                    lcs.nodes[lcs.index[i].tail].match.col + pos - 2 * i,
                        .link = lcs.nodes[lcs.index[i].tail].idx,
                    }));

                prev = lcs.index[i].tail;
            } // if

            graph.dom_nodes[array_length(graph.dom_nodes) - 1].match.length += 1;
        } // if
    } // for

    // add empty match sink if not already handled because of unique matches
    if (graph.dom_nodes[array_length(graph.dom_nodes) - 1].link != sink.idx)
    {
        ARRAY_APPEND(allocator, graph.dom_nodes,
            ((GVA_Dom_Node)
            {
                .match = {sink.match.row + sink.match.length, sink.match.col + sink.match.length},
                .distance = distance,
                .link = sink.idx,
            }));
    } // if

    // update source and sink according to the supremal
    if (graph.dom_nodes != NULL)
    {
        graph.nodes[source.idx].match.row += graph.dom_nodes[0].match.length;
        graph.nodes[source.idx].match.col += graph.dom_nodes[0].match.length;
        graph.nodes[source.idx].match.length -= graph.dom_nodes[0].match.length;

        graph.dom_nodes[0].match.row += graph.dom_nodes[0].match.length;
        graph.dom_nodes[0].match.col += graph.dom_nodes[0].match.length;
        graph.dom_nodes[0].match.length = 0;

        graph.nodes[sink.idx].match.length -= graph.dom_nodes[array_length(graph.dom_nodes) - 1].match.length;

        graph.dom_nodes[array_length(graph.dom_nodes) - 1].match.length = 0;
    } // if

    graph.source = source.idx;

    lcs.index = allocator.allocate(allocator.context, lcs.index, lcs.length * sizeof(*lcs.index), 0);
    lcs.nodes = ARRAY_DESTROY(allocator, lcs.nodes);
    return graph;
} // gva_lcs_graph_init


inline void
gva_lcs_graph_destroy(GVA_Allocator const allocator,
    GVA_LCS_Graph self, bool const observed)
{
    self.nodes = ARRAY_DESTROY(allocator, self.nodes);
    self.edges = ARRAY_DESTROY(allocator, self.edges);
    self.dom_nodes = ARRAY_DESTROY(allocator, self.dom_nodes);
    if (observed)
    {
        gva_string_destroy(allocator, self.observed);
    } // if
} // gva_lcs_graph_destroy


inline GVA_LCS_Graph
gva_lcs_graph_from_allele(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const n, GVA_Variant const variants[static restrict n])
{
    GVA_LCS_Graph graph = {.observed = gva_patch(allocator, len_ref, reference, n, variants)};
    local_supremal(allocator, len_ref, reference, graph.observed.len, graph.observed.str, 0, 0, &graph);
    if (graph.nodes == NULL)
    {
        graph = gva_lcs_graph_init(allocator, len_ref, reference, graph.observed.len, graph.observed.str, 0);
    } // if
    return graph;
} // gva_lcs_graph_from_allele


GVA_LCS_Graph
gva_lcs_graph_from_variants(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const n, GVA_Variant const variants[static restrict n])
{
    GVA_Variant variant = {.start = variants[0].start, .end = variants[n - 1].end};
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
            return (GVA_LCS_Graph) {NULL};
        } // if

        memcpy(observed, reference + start, variant.start - start);
        memcpy(observed + variant.start - start, variant.sequence.str, variant.sequence.len);
        memcpy(observed + variant.start - start + variant.sequence.len, reference + variant.end, end - variant.end);

        GVA_LCS_Graph graph = gva_lcs_graph_init(allocator, end - start, reference + start, len, observed, start);
        GVA_Variant const supremal = gva_lcs_graph_supremal(graph);

        if ((supremal.start > start || supremal.start == 0) &&
            (supremal.end   < end   || supremal.end   == len_ref))
        {
            // observed is now owned by graph
            gva_string_destroy(allocator, variant.sequence);
            return graph;
        } // if

        gva_lcs_graph_destroy(allocator, graph, false);

        old_len = len;
        offset *= 2;  // OVERFLOW
    } // while
} // gva_lcs_graph_from_variants


inline size_t
gva_lcs_graph_distance(GVA_LCS_Graph const self)
{
    return self.dom_nodes[array_length(self.dom_nodes) - 1].distance;
} // gva_lcs_graph_distance


inline GVA_Variant
gva_lcs_graph_local_supremal(GVA_LCS_Graph const self,
    size_t const start, size_t const end)
{
    GVA_Variant variant;
    gva_edges(self.observed.str,
              self.dom_nodes[start].match, self.dom_nodes[end].match,
              start == 0, end == array_length(self.dom_nodes) - 1,
              &variant);
    return variant;
} // gva_lcs_graph_local_supremal


inline GVA_Variant
gva_lcs_graph_supremal(GVA_LCS_Graph const self)
{
    return gva_lcs_graph_local_supremal(self, 0, array_length(self.dom_nodes) - 1);
} // gva_lcs_graph_supremal


inline size_t
gva_edges(char const observed[static restrict 1],
    GVA_Match const head, GVA_Match const tail,
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
