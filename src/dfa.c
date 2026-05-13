#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint8_t, SIZE_MAX
#include <string.h>     // memset

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_edges
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "../include/variant.h"     // GVA_Variant
#include "align.h"              // LCS_Alignment, lcs_align
#include "array.h"              // ARRAY_*, array_*
#include "common.h"             // MAX, MIN
#include "dfa.h"                // DFAs, dfa_*
#include "hash_table.h"         // HASH_TABLE_*, hash_table_*


// FIXME: DEBUG
#include <stdio.h>


static uint8_t const DFA_DELETION  = 0x1;
static uint8_t const DFA_INSERTION = 0x2;


static inline uint8_t
get(uint8_t const dfa[static 1], size_t const idx)
{
    return dfa[idx / 4] >> (2 * (idx % 4));
} // get


static inline void
set(uint8_t dfa[static 1], size_t const idx, uint8_t const value)
{
    dfa[idx / 4] |= value << (2 * (idx % 4));
} // set


static inline void
edge(uint8_t dfa[static 1], size_t const width,
    size_t const start, size_t const end,
    size_t const len, size_t const offset)
{
    for (size_t i = start; i <= end; ++i)
    {
        for (size_t j = 0; j <= len; ++j)
        {
            set(dfa, i * width + offset + j,
                (i < end) * DFA_DELETION | (j < len) * DFA_INSERTION);
        } // for
    } // for
} // edge


uint8_t*
dfa_concat(GVA_Allocator const allocator, size_t const n,
    GVA_Variant const variants[static restrict n],
    uint8_t* const dfas[static restrict n])
{
    if (n == 0)
    {
        return NULL;
    } // if
    if (n == 1)
    {
        return dfas[0];
    } // if

    size_t const start = variants[0].start;
    size_t const end = variants[n - 1].end;

    size_t width = 1;
    for (size_t i = 0; i < n; ++i)
    {
        if (i > 0)
        {
            width += variants[i].start - variants[i - 1].end;
        } // if
        width += variants[i].sequence.len;
    } // for
    size_t const size = (end - start + 1) * width;
    size_t const bytes = (size + 3) / 4;

    uint8_t* const restrict dfa = allocator.allocate(allocator.context, NULL, 0, bytes);
    if (dfa == NULL)
    {
        return NULL;  // OOM
    } // if

    memset(dfa, 0, bytes);

    size_t offset_col = 0;
    for (size_t i = 0; i < n; ++i)
    {
        size_t const offset_row = variants[i].start - start;
        if (i > 0)
        {
            offset_col += variants[i - 1].sequence.len + variants[i].start - variants[i - 1].end;
        } // if
        for (size_t j = 0; j <= variants[i].end - variants[i].start; ++j)
        {
            for (size_t k = 0; k <= variants[i].sequence.len; ++k)
            {
                fprintf(stderr, "cpy %zu, %zu -> %zu, %zu\n", j, k, j + offset_row, k + offset_col);
                set(dfa, (j + offset_row) * width + k + offset_col,
                    get(dfas[i], j * (variants[i].sequence.len + 1) + k));
            } // for
        } // for
    } // for

    return dfa;
} // dfa_concat


uint8_t*
dfa_from_alignment(GVA_Allocator const allocator, uint8_t* restrict dfas,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs])
{
    LCS_Alignment lcs = lcs_align(allocator, len_ref, reference, len_obs, observed, 0);
    size_t const distance = len_ref + len_obs - 2 * lcs.length;

    if (distance == 0)
    {
        return dfas;
    } // if

    size_t const width = len_obs + 1;
    size_t const size = (len_ref + 1) * width;
    size_t const bytes = (size + 3) / 4;

    dfas = array_ensure(allocator, dfas, 1, bytes);
    if (dfas == NULL)
    {
        lcs.index = allocator.allocate(allocator.context, lcs.index, lcs.length * sizeof(*lcs.index), 0);
        lcs.nodes = ARRAY_DESTROY(allocator, lcs.nodes);
        return NULL;  // OOM
    } // if

    size_t const start = array_length(dfas);
    memset(dfas + start, 0, bytes);
    array_header(dfas)->length += bytes;

    if (lcs.nodes == NULL)
    {
        edge(dfas + start, width, 0, len_ref, len_obs, 0);
        return dfas;
    } // if

    gva_uint tail_idx = lcs.index[lcs.length - 1].tail;
    LCS_Node sink = lcs.nodes[tail_idx];
    if (sink.match.row + sink.match.length == len_ref &&
        sink.match.col + sink.match.length == len_obs)
    {
        lcs.nodes[tail_idx].moved = true;
        sink = lcs.nodes[tail_idx];
    } // if
    else
    {
        sink = (LCS_Node) {.match = {len_ref, len_obs}};
        tail_idx = GVA_NULL;
    } // else
    for (gva_uint i = lcs.index[lcs.length - 1].head; i != tail_idx; i = lcs.nodes[i].next)
    {
        edge(dfas + start, width,
            lcs.nodes[i].match.row + lcs.nodes[i].match.length, sink.match.row + sink.match.length,
            sink.match.col + sink.match.length - (lcs.nodes[i].match.col + lcs.nodes[i].match.length),
            lcs.nodes[i].match.col + lcs.nodes[i].match.length);
        lcs.nodes[i].moved = true;
    } // for

    for (gva_uint i = lcs.length - 1; i >= 1; --i)
    {
        gva_uint next = GVA_NULL;
        for (gva_uint j = lcs.index[i].head; j != GVA_NULL; j = next)
        {
            LCS_Node* const restrict tail = &lcs.nodes[j];
            next = tail->next;
            if (!tail->moved)
            {
                continue;
            } // if

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

                edge(dfas + start, width,
                    head->match.row + head->match.length, tail->match.row + tail->match.length - 1,
                    tail->match.col + tail->match.length - (head->match.col + head->match.length) - 1,
                    head->match.col + head->match.length);
                head->moved = true;

                here = k;
            } // for

            if (tail->match.length > 1)
            {
                tail->match.length -= 1;
                if (here != GVA_NULL)
                {
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
    if (source.match.row == 0 && source.match.col == 0)
    {
        head_idx = source.next;
    } // if
    else
    {
        source = (LCS_Node) {0};
    } // else
    for (gva_uint i = head_idx; i != GVA_NULL; i = lcs.nodes[i].next)
    {
        if (!lcs.nodes[i].moved)
        {
            continue;
        } // if

        edge(dfas + start, width, source.match.row, lcs.nodes[i].match.row,
            lcs.nodes[i].match.col - source.match.col, 0);
    } // for

    lcs.index = allocator.allocate(allocator.context, lcs.index, lcs.length * sizeof(*lcs.index), 0);
    lcs.nodes = ARRAY_DESTROY(allocator, lcs.nodes);

    return dfas;
} // dfa_from_alignment


DFAs
dfas_from_lcs_graph(GVA_Allocator const allocator, GVA_LCS_Graph const graph)
{
    size_t bytes = 0;
    for (size_t i = 0; i < array_length(graph.dom_nodes) - 1; ++i)
    {
        GVA_Variant const local = gva_lcs_graph_local_supremal(graph, i, i + 1);
        bytes += ((local.end - local.start + 1) * (local.sequence.len + 1) + 3) / 4;
    } // for

    DFAs dfas = 
    {
        .len = array_length(graph.dom_nodes),
        .starts = allocator.allocate(allocator.context, NULL, 0, sizeof(*dfas.starts) * array_length(graph.dom_nodes)),
        .data = allocator.allocate(allocator.context, NULL, 0, bytes),
    };

    gva_uint* restrict next = allocator.allocate(allocator.context, NULL, 0, sizeof(*next) * array_length(graph.nodes));

    if (dfas.starts == NULL || dfas.data == NULL || next == NULL)
    {
        next = allocator.allocate(allocator.context, next, sizeof(*next) * array_length(graph.nodes), 0);
        dfas.data = allocator.allocate(allocator.context, dfas.data, bytes, 0);
        dfas.starts = allocator.allocate(allocator.context, dfas.starts, sizeof(*dfas.starts) * array_length(graph.dom_nodes), 0);
        return (DFAs) {0};  // OOM
    } // if

    memset(dfas.data, 0, bytes);
    for (size_t i = 0; i < array_length(graph.nodes); ++i)
    {
        next[i] = GVA_NULL;
    } // for

    fprintf(stderr, "%zu bytes\n", bytes);

    size_t start = 0;
    for (size_t i = 0; i < array_length(graph.dom_nodes) - 1; ++i)
    {
        GVA_Variant const local = gva_lcs_graph_local_supremal(graph, i, i + 1);
        size_t const width = local.sequence.len + 1;

        fprintf(stderr, "DFA %u x %zu @ %zu\n", local.end - local.start + 1, width, start);

        gva_uint head = graph.dom_nodes[i].link;
        gva_uint tail = head;
        while (head != GVA_NULL)
        {
            fprintf(stderr, "POP %u\n", head);

            if (head == graph.dom_nodes[i + 1].link)
            {
                head = next[head];
                continue;
            } // if

            for (gva_uint j = graph.nodes[head].edges; j != GVA_NULL; j = graph.edges[j].next)
            {
                GVA_Variant variant = {0};
                size_t const count = gva_edges(graph.observed.str,
                    graph.nodes[head].match, graph.nodes[graph.edges[j].tail].match,
                    head == graph.source, graph.nodes[graph.edges[j].tail].edges == GVA_NULL,
                    &variant);
                for (size_t k = 0; k < count; ++k)
                {
                    edge(dfas.data + start, width, variant.start + k - local.start, variant.end + k - local.start,
                        variant.sequence.len, variant.sequence.str - local.sequence.str + k);
                } // for

                if (next[graph.edges[j].tail] == GVA_NULL)
                {
                    fprintf(stderr, "  PUSH %u\n", graph.edges[j].tail);
                    next[tail] = graph.edges[j].tail;
                    tail = graph.edges[j].tail;
                } // if
            } // for

            head = next[head];
        } // while
        dfas.starts[i] = start;
        start += ((local.end - local.start + 1) * width + 3) / 4;
    } // for
    dfas.starts[array_length(graph.dom_nodes) - 1] = bytes;

    next = allocator.allocate(allocator.context, next, sizeof(*next) * array_length(graph.nodes), 0);

    return dfas;
} // dfas_from_lcs_graph


typedef struct
{
    HASH_TABLE_KEY;
    gva_uint included;
    size_t   next;
} Entry;


typedef struct
{
    Entry*        entries;
    GVA_Allocator allocator;
    size_t        head;
    size_t        tail;
} Queue;


static inline Queue
queue_init(GVA_Allocator const allocator, size_t const capacity)
{
    return (Queue)
    {
        .entries = hash_table_init(allocator, capacity, sizeof(Entry)),
        .allocator = allocator,
        .head = SIZE_MAX,
        .tail = SIZE_MAX,
    };
} // queue_init


static inline void
queue_destroy(Queue self[static 1])
{
    self->entries = HASH_TABLE_DESTROY(self->allocator, self->entries);
    self->head = SIZE_MAX;
    self->tail = SIZE_MAX;
} // queue_destroy


static inline bool
queue_empty(Queue const self[static 1])
{
    return self->head == SIZE_MAX;
} // queue_empty


static inline size_t
queue_pop(Queue self[static 1])
{
    size_t const hash_idx = HASH_TABLE_INDEX(self->entries, self->head);
    size_t const idx = self->entries[hash_idx].gva_key;
    self->head = self->entries[hash_idx].next;
    if (self->head == SIZE_MAX)
    {
        self->tail = SIZE_MAX;
    } // if
    return idx;
} // queue_pop


static inline void
queue_push_back(Queue self[static 1], size_t const idx, size_t const included)
{
    size_t const hash_idx = HASH_TABLE_INDEX(self->entries, idx);
    if (self->entries[hash_idx].gva_key == idx)
    {
        return;
    } // if

    HASH_TABLE_SET(self->allocator, self->entries, idx,
        ((Entry)
        {
            .gva_key = idx,
            .included = included,
            .next = SIZE_MAX,
        }));

    if (self->tail == SIZE_MAX)
    {
        self->head = idx;
    } // if
    else
    {
        self->entries[HASH_TABLE_INDEX(self->entries, self->tail)].next = idx;
    } // else
    self->tail = idx;
} // queue_push_back


static inline void
queue_push_front(Queue self[static 1], size_t const idx, size_t const included)
{
    size_t const hash_idx = HASH_TABLE_INDEX(self->entries, idx);
    if (self->entries[hash_idx].gva_key == idx)
    {
        return;
    } // if

    HASH_TABLE_SET(self->allocator, self->entries, idx,
        ((Entry)
        {
            .gva_key = idx,
            .included = included,
            .next = self->head,
        }));

    if (self->head == SIZE_MAX)
    {
        self->tail = idx;
    } // if
    self->head = idx;
} // queue_push_front


size_t
dfa_max_overlap(GVA_Allocator const allocator,
    size_t const len, char const reference[static restrict len],
    GVA_Variant const lhs_variant, uint8_t const lhs_dfa[static restrict 1],
    GVA_Variant const rhs_variant, uint8_t const rhs_dfa[static restrict 1],
    size_t const limit)
{
    static size_t const INITIAL_SIZE = 4096;

    size_t const lhs_width = lhs_variant.sequence.len + 1;
    size_t const rhs_width = rhs_variant.sequence.len + 1;
    size_t const lhs_size = (lhs_variant.end - lhs_variant.start + 1) * lhs_width;
    size_t const rhs_size = (rhs_variant.end - rhs_variant.start + 1) * rhs_width;

    Queue queue = queue_init(allocator, INITIAL_SIZE);
    if (queue.entries == NULL)
    {
        return 0;  // OOM
    } // if

    //size_t count = 0;
    queue_push_front(&queue, 0, 0);
    while (!queue_empty(&queue))
    {
        size_t const idx = queue_pop(&queue);
        if (idx == lhs_size * rhs_size - 1)  // OVERFLOW
        {
            break;
        } // if

        size_t const lhs_idx = idx / rhs_size;
        size_t const rhs_idx = idx % rhs_size;
        size_t const lhs_ref = lhs_idx / lhs_width + lhs_variant.start;
        size_t const rhs_ref = rhs_idx / rhs_width + rhs_variant.start;
        size_t const included = queue.entries[HASH_TABLE_INDEX(queue.entries, idx)].included;

        //count += 1;
        // fprintf(stderr, "{%zu, %zu} (%zu): %zu\n", lhs_idx, rhs_idx, idx, included);

        if ((lhs_idx == 0 && rhs_ref < lhs_ref) || lhs_idx == lhs_size - 1)
        {
            uint8_t const rhs_value = get(rhs_dfa, rhs_idx);
            bool const rhs_match = rhs_idx % rhs_width < rhs_width - 1 &&
                reference[rhs_variant.start + rhs_idx / rhs_width] == rhs_variant.sequence.str[rhs_idx % rhs_width];

            if (rhs_match)
            {
                // ..
                // fprintf(stderr, "    push (. vs %zu MU) +0 +0 {%zu, %zu} (%zu)\n", rhs_ref, lhs_idx, rhs_idx + rhs_width + 1, lhs_idx * rhs_size + rhs_idx + rhs_width + 1);
                queue_push_front(&queue, lhs_idx * rhs_size + rhs_idx + rhs_width + 1, included);
                continue;
            } // if

            bool const rhs_deletion = rhs_value & DFA_DELETION;
            bool const rhs_insertion = rhs_value & DFA_INSERTION;

            if (rhs_deletion)
            {
                // .D
                // fprintf(stderr, "    push (. vs %zu DELTA) +0 +1 {%zu, %zu} (%zu)\n", rhs_ref, lhs_idx, rhs_idx + rhs_width, lhs_idx * rhs_size + rhs_idx + rhs_width);
                queue_push_back(&queue, lhs_idx * rhs_size + rhs_idx + rhs_width, included);
            } // if
            if (rhs_insertion)
            {
                // .I
                // fprintf(stderr, "    push (. vs %zu %c) +0 +1 {%zu, %zu} (%zu)\n", rhs_ref, rhs_variant.sequence.str[rhs_idx % rhs_width], lhs_idx, rhs_idx + 1, lhs_idx * rhs_size + rhs_idx + 1);
                queue_push_back(&queue, lhs_idx * rhs_size + rhs_idx + 1, included);
            } // if
            continue;
        } // if

        if ((rhs_idx == 0 && lhs_ref < rhs_ref) || rhs_idx == rhs_size - 1)
        {
            uint8_t const lhs_value = get(lhs_dfa, lhs_idx);
            bool const lhs_match = lhs_idx % lhs_width < lhs_width - 1 &&
                reference[lhs_variant.start + lhs_idx / lhs_width] == lhs_variant.sequence.str[lhs_idx % lhs_width];

            if (lhs_match)
            {
                // ..
                // fprintf(stderr, "    push (%zu MU vs .) +0 +0 {%zu, %zu} (%zu)\n", lhs_ref, lhs_idx + lhs_width + 1, rhs_idx, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx);
                queue_push_front(&queue, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx, included);
                continue;
            } // if

            bool const lhs_deletion = lhs_value & DFA_DELETION;
            bool const lhs_insertion = lhs_value & DFA_INSERTION;

            if (lhs_deletion)
            {
                // D.
                // fprintf(stderr, "    push (%zu DELTA vs .) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, lhs_idx + lhs_width, rhs_idx, (lhs_idx + lhs_width) * rhs_size + rhs_idx);
                queue_push_back(&queue, (lhs_idx + lhs_width) * rhs_size + rhs_idx, included);
            } // if
            if (lhs_insertion)
            {
                // I.
                // fprintf(stderr, "    push (%zu %c vs .) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, lhs_variant.sequence.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx, (lhs_idx + 1) * rhs_size + rhs_idx);
                queue_push_back(&queue, (lhs_idx + 1) * rhs_size + rhs_idx, included);
            } // if
            continue;
        } // if

        uint8_t const lhs_value = get(lhs_dfa, lhs_idx);
        bool const lhs_match = lhs_idx % lhs_width < lhs_width - 1 &&
            reference[lhs_variant.start + lhs_idx / lhs_width] == lhs_variant.sequence.str[lhs_idx % lhs_width];
        uint8_t const rhs_value = get(rhs_dfa, rhs_idx);
        bool const rhs_match = rhs_idx % rhs_width < rhs_width - 1 &&
            reference[rhs_variant.start + rhs_idx / rhs_width] == rhs_variant.sequence.str[rhs_idx % rhs_width];

        if (lhs_match && rhs_match)
        {
            // ..
            // fprintf(stderr, "    push (%zu MU) {%zu, %zu} +0 +0 (%zu)\n", lhs_ref, lhs_idx + lhs_width + 1, rhs_idx + rhs_width + 1, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx + rhs_width + 1);
            queue_push_front(&queue, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx + rhs_width + 1, included);
            continue;
        } // if

        bool const lhs_deletion = lhs_value & DFA_DELETION;
        bool const lhs_insertion = lhs_value & DFA_INSERTION;
        bool const rhs_deletion = rhs_value & DFA_DELETION;
        bool const rhs_insertion = rhs_value & DFA_INSERTION;

        if (lhs_deletion && rhs_deletion)
        {
            // DD
            // fprintf(stderr, "    push (%zu DELTA) {%zu, %zu} +1 +0 (%zu)\n", lhs_ref, lhs_idx + lhs_width, rhs_idx + rhs_width, (lhs_idx + lhs_width) * rhs_size + rhs_idx + rhs_width);
            if (limit > 0 && included + 1 >= limit)
            {
                break;
            } // if
            queue_push_front(&queue, (lhs_idx + lhs_width) * rhs_size + rhs_idx + rhs_width, included + 1);
        } // if

        if (lhs_insertion && rhs_insertion && lhs_variant.sequence.str[lhs_idx % lhs_width] == rhs_variant.sequence.str[rhs_idx % rhs_width])
        {
            // II
            // fprintf(stderr, "    push (%zu %c) {%zu, %zu} +1 +0 (%zu)\n", lhs_ref, lhs_variant.sequence.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx + 1, (lhs_idx + 1) * rhs_size + rhs_idx + 1);
            if (limit > 0 && included + 1 >= limit)
            {
                break;
            } // if
            queue_push_front(&queue, (lhs_idx + 1) * rhs_size + rhs_idx + 1, included + 1);
        } // if

        if (lhs_deletion && rhs_match)
        {
            // D.
            // fprintf(stderr, "    push (%zu DELTA vs %zu MU) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, rhs_ref, lhs_idx + lhs_width, rhs_idx + rhs_width + 1, (lhs_idx + lhs_width) * rhs_size + rhs_idx + rhs_width + 1);
            queue_push_back(&queue, (lhs_idx + lhs_width) * rhs_size + rhs_idx + rhs_width + 1, included);
        } // if
        else if (lhs_match && rhs_deletion)
        {
            // .D
            // fprintf(stderr, "    push (%zu MU vs %zu DELTA) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, rhs_ref, lhs_idx + lhs_width + 1, rhs_idx + rhs_width, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx + rhs_width);
            queue_push_back(&queue, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx + rhs_width, included);
        } // if

        if (lhs_insertion && (rhs_deletion || rhs_match || (rhs_insertion &&
            lhs_variant.sequence.str[lhs_idx % lhs_width] != rhs_variant.sequence.str[rhs_idx % rhs_width])))
        {
            // I.
            // fprintf(stderr, "    push (%zu %c vs .) {%zu, %zu} +0 +1 (%zu)\n", lhs_ref, lhs_variant.sequence.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx, (lhs_idx + 1) * rhs_size + rhs_idx);
            queue_push_back(&queue, (lhs_idx + 1) * rhs_size + rhs_idx, included);
        } // if

        if (rhs_insertion && (lhs_deletion || lhs_match || (lhs_insertion &&
            lhs_variant.sequence.str[lhs_idx % lhs_width] != rhs_variant.sequence.str[rhs_idx % rhs_width])))
        {
            // .I
            // fprintf(stderr, "    push (. vs %zu %c) {%zu, %zu} +0 +1 (%zu)\n", rhs_ref, rhs_variant.sequence.str[rhs_idx % rhs_width], lhs_idx, rhs_idx + 1, lhs_idx * rhs_size + rhs_idx + 1);
            queue_push_back(&queue, lhs_idx * rhs_size + rhs_idx + 1, included);
        } // if
    } // while

    //fprintf(stderr, "count: %zu\n", count);

    size_t included = limit;
    size_t const idx = HASH_TABLE_INDEX(queue.entries, lhs_size * rhs_size - 1);
    if (queue.entries[idx].gva_key == lhs_size * rhs_size - 1)
    {
        included = queue.entries[idx].included;
    } // if

    queue_destroy(&queue);

    return included;
} // dfa_max_overlap


bool
dfa_disjoint(GVA_Variant const lhs_variant, uint8_t const lhs_dfa[static restrict 1],
    GVA_Variant const rhs_variant, uint8_t const rhs_dfa[static restrict 1])
{
    size_t const lhs_width = lhs_variant.sequence.len + 1;
    size_t const rhs_width = rhs_variant.sequence.len + 1;

    size_t const start = MAX(lhs_variant.start, rhs_variant.start);
    size_t const end = MIN(lhs_variant.end + 1, rhs_variant.end + 1);

    for (size_t i = start; i < end; ++i)
    {
        for (size_t j = 0; j < lhs_width; ++j)
        {
            uint8_t const lhs_value = get(lhs_dfa, (i - lhs_variant.start) * lhs_width + j);

            for (size_t k = 0; k < rhs_width; ++k)
            {
                uint8_t const rhs_value = get(rhs_dfa, (i - rhs_variant.start) * rhs_width + k);

                if ((lhs_value & DFA_DELETION && rhs_value & DFA_DELETION) ||
                    (lhs_value & DFA_INSERTION && rhs_value & DFA_INSERTION &&
                    lhs_variant.sequence.str[j] == rhs_variant.sequence.str[k]))
                {
                    return false;
                } // if
            } // for
        } // for
    } // for

    return true;
} // dfa_disjoint


// FIXME: DEBUG
void
dfa_dot(FILE* restrict const stream,
    size_t const len, char const reference[static restrict len],
    GVA_Variant const variant, uint8_t const dfa[static restrict 1])
{
    size_t const width = variant.sequence.len + 1;
    size_t const size = (variant.end - variant.start + 1) * width;

    fprintf(stream, "strict digraph{\nrankdir=LR\nnode[fixedsize=true,shape=circle,width=1]\ni[label=\"\",shape=none,width=0]\ni->0\n");

    size_t end = 0;
    size_t j = 0;
    for (size_t i = 0; i <= variant.end - variant.start; ++i)
    {
        size_t start = width;
        size_t next_end = end;
        uint8_t value = 0;
        while (j <= end || value & DFA_INSERTION)
        {
            value = get(dfa, i * width + j);

            if (value & DFA_DELETION)
            {
                fprintf(stream, "%zu->%zu[label=\"%zu &Delta;\"]\n", i * width + j, i * width + j + width, variant.start + i);
                start = MIN(start, j);
            } // if

            if (value & DFA_INSERTION)
            {
                fprintf(stream, "%zu->%zu[label=\"%zu %c\"]\n", i * width + j, i * width + j + 1, variant.start + i, variant.sequence.str[j]);
                next_end = MAX(next_end, j);
            } // if

            if (i < variant.end - variant.start &&
                reference[variant.start + i] == variant.sequence.str[j])
            {
                fprintf(stream, "%zu->%zu[label=\"%zu &Mu;\"]\n", i * width + j, i * width + j + width + 1, variant.start + i);
                start = MIN(start, j + 1);
                next_end = MAX(next_end, j + 1);
            } // if

            j += 1;
        } // while
        end = next_end;
        j = start;
    } // for

    fprintf(stream, "%zu[peripheries=2]\n}\n", size - 1);
} // dfa_dot


void
dfa_svg(FILE* restrict const stream,
    size_t const len, char const reference[static restrict len],
    GVA_Variant const variant, uint8_t const dfa[static restrict 1])
{
    size_t const width = variant.sequence.len + 1;
    size_t const size = (variant.end - variant.start + 1) * width;

    fprintf(stream, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 %zu %u\">\n", width * 100, (variant.end - variant.start + 1) * 100);
    fprintf(stream, "<defs>\n<marker id=\"a\" markerHeight=\"8\" markerWidth=\"10\" orient=\"auto\" refX=\"2\" refY=\"4\">\n<path fill=\"context-stroke\" d=\"m0 0 2 4-2 4 10-4Z\"/>\n</marker>\n</defs>\n");
    fprintf(stream, "<style>\n:root{--c:#1f2328;}\ntext{dominant-baseline:middle;fill:var(--c);font-size:smaller;text-anchor:middle;}\ncircle{fill:none;stroke:var(--c);}\nline{marker-end:url(#a);stroke:var(--c);}\n</style>\n");
    fprintf(stream, "<g font-size=\"larger\">\n");
    fprintf(stream, "<text x=\"7\" y=\"12\">%u</text>\n", variant.start);
    for (size_t i = 0; i < width - 1; ++i)
    {
        fprintf(stream, "<text x=\"%zu\" y=\"12\">%c</text>\n", (i + 1) * 100, variant.sequence.str[i]);
    } // for
    fprintf(stream, "</g>\n");

    size_t end = 0;
    size_t j = 0;
    for (size_t i = 0; i <= variant.end - variant.start; ++i)
    {
        size_t start = width;
        size_t next_end = end;
        uint8_t value = 0;
        while (j <= end || value & DFA_INSERTION)
        {
            value = get(dfa, i * width + j);

            bool edge = false;
            if (value & DFA_DELETION)
            {
                fprintf(stream, "<line x1=\"%zu\" y1=\"%zu\" x2=\"%zu\" y2=\"%zu\" />\n", j * 100 + 50, i * 100 + 50 + 25, j * 100 + 50, i * 100 + 100 + 16);
                start = MIN(start, j);
                edge = true;
            } // if

            if (value & DFA_INSERTION)
            {
                fprintf(stream, "<line x1=\"%zu\" y1=\"%zu\" x2=\"%zu\" y2=\"%zu\" />\n", j * 100 + 50 + 25, i * 100 + 50, j * 100 + 100 + 16, i * 100 + 50);
                next_end = MAX(next_end, j);
                edge = true;
            } // if

            if (i < variant.end - variant.start &&
                reference[variant.start + i] == variant.sequence.str[j])
            {
                fprintf(stream, "<line x1=\"%zu\" y1=\"%zu\" x2=\"%zu\" y2=\"%zu\" />\n", j * 100 + 50 + 18, i * 100 + 50 + 18, j * 100 + 100 + 50 - 24, i * 100 + 100 + 50 - 24);
                start = MIN(start, j + 1);
                next_end = MAX(next_end, j + 1);
                edge = true;
            } // if

            if (edge)
            {
                fprintf(stream, "<circle cx=\"%zu\" cy=\"%zu\" r=\"25\" />\n", j * 100 + 50, i * 100 + 50);
                fprintf(stream, "<text x=\"%zu\" y=\"%zu\">%zu</text>\n", j * 100 + 50, i * 100 + 50, i * width + j);
            } // if

            j += 1;
        } // while
        end = next_end;
        j = start;
    } // for

    fprintf(stream, "<line x1=\"0\" y1=\"50\" x2=\"16\" y2=\"50\" />\n");
    fprintf(stream, "<circle cx=\"%zu\" cy=\"%zu\" r=\"25\" />\n", ((size - 1) % width) * 100 + 50, ((size - 1) / width) * 100 + 50);
    fprintf(stream, "<circle cx=\"%zu\" cy=\"%zu\" r=\"22\" />\n", ((size - 1) % width) * 100 + 50, ((size - 1) / width) * 100 + 50);
                fprintf(stream, "<text x=\"%zu\" y=\"%zu\">%zu</text>\n", ((size - 1) % width) * 100 + 50, ((size - 1) / width) * 100 + 50, size - 1);
    fprintf(stream, "</svg>\n");
} // dfa_svg
