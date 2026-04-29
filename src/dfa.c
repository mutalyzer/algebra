#include <stdbool.h>    // bool, true
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint8_t, SIZE_MAX
#include <string.h>     // memset

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_edges
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "../include/variant.h"     // GVA_Variant
#include "align.h"              // LCS_Alignment, lcs_align
#include "array.h"              // ARRAY_*, array_*
#include "dfa.h"                // dfa*
#include "hash_table.h"         // HASH_TABLE_*, hash_table_*


#include <stdio.h>      // DEBUG


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
dfa_from_alignment(GVA_Allocator const allocator, uint8_t* dfas,
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
    size_t const len = (size + 3) / 4;

    dfas = array_ensure(allocator, dfas, 1, len);
    if (dfas == NULL)
    {
        lcs.index = allocator.allocate(allocator.context, lcs.index, lcs.length * sizeof(*lcs.index), 0);
        lcs.nodes = ARRAY_DESTROY(allocator, lcs.nodes);
        return NULL;  // OOM
    } // if

    size_t const start = array_length(dfas);
    memset(dfas + start, 0, len);
    array_header(dfas)->length += len;

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


uint8_t*
dfa_from_lcs_graph(GVA_Allocator const allocator, uint8_t* dfas,
    GVA_LCS_Graph const graph)
{
    GVA_Variant const supremal = gva_lcs_graph_supremal(graph);
    size_t const width = supremal.sequence.len + 1;
    size_t const size = (supremal.end - supremal.start + 1) * width;
    size_t const len = (size + 3) / 4;

    dfas = array_ensure(allocator, dfas, 1, len);
    if (dfas == NULL)
    {
        return NULL;  // OOM
    } // if

    size_t const start = array_length(dfas);
    memset(dfas + start, 0, len);
    array_header(dfas)->length += len;

    for (size_t i = 0; i < array_length(graph.nodes); ++i)
    {
        for (gva_uint j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
        {
            GVA_Variant variant = {0};
            size_t const count = gva_edges(graph.observed.str,
                graph.nodes[i].match, graph.nodes[graph.edges[j].tail].match,
                i == graph.source, graph.nodes[graph.edges[j].tail].edges == GVA_NULL,
                &variant);
            for (size_t k = 0; k < count; ++k)
            {
                edge(dfas + start, width, variant.start + k - supremal.start, variant.end + k - supremal.start,
                    variant.sequence.len, variant.sequence.str - supremal.sequence.str + k);
            } // for
        } // for
    } // for

    return dfas;
} // dfa_from_lcs_graph


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
    GVA_Variant const lhs_supremal, uint8_t const lhs_dfa[static restrict 1],
    GVA_Variant const rhs_supremal, uint8_t const rhs_dfa[static restrict 1])
{
    static size_t const INITIAL_SIZE = 4096;

    size_t const lhs_width = lhs_supremal.sequence.len + 1;
    size_t const rhs_width = rhs_supremal.sequence.len + 1;
    size_t const lhs_size = (lhs_supremal.end - lhs_supremal.start + 1) * lhs_width;
    size_t const rhs_size = (rhs_supremal.end - rhs_supremal.start + 1) * rhs_width;

    Queue queue = queue_init(allocator, INITIAL_SIZE);
    if (queue.entries == NULL)
    {
        return 0;  // OOM
    } // if

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
        size_t const lhs_ref = lhs_idx / lhs_width + lhs_supremal.start;
        size_t const rhs_ref = rhs_idx / rhs_width + rhs_supremal.start;
        size_t const included = queue.entries[HASH_TABLE_INDEX(queue.entries, idx)].included;

        //fprintf(stderr, "{%zu, %zu} (%zu): %zu\n", lhs_idx, rhs_idx, idx, included);

        if ((lhs_idx == 0 && rhs_ref < lhs_ref) || lhs_idx == lhs_size - 1)
        {
            uint8_t const rhs_value = get(rhs_dfa, rhs_idx);
            bool const rhs_deletion = rhs_value & DFA_DELETION;
            bool const rhs_insertion = rhs_value & DFA_INSERTION;
            bool const rhs_match = rhs_idx % rhs_width < rhs_width - 1 &&
                reference[rhs_supremal.start + rhs_idx / rhs_width] == rhs_supremal.sequence.str[rhs_idx % rhs_width];

            if (rhs_match)
            {
                // ..
                // fprintf(stderr, "    push (. vs %zu MU) +0 +0 {%zu, %zu} (%zu)\n", rhs_ref, lhs_idx, rhs_idx + rhs_width + 1, lhs_idx * rhs_size + rhs_idx + rhs_width + 1);
                queue_push_front(&queue, lhs_idx * rhs_size + rhs_idx + rhs_width + 1, included);
            } // if

            if (rhs_deletion)
            {
                // .D
                // fprintf(stderr, "    push (. vs %zu DELTA) +0 +1 {%zu, %zu} (%zu)\n", rhs_ref, lhs_idx, rhs_idx + rhs_width, lhs_idx * rhs_size + rhs_idx + rhs_width);
                queue_push_back(&queue, lhs_idx * rhs_size + rhs_idx + rhs_width, included);
            } // if

            if (rhs_insertion)
            {
                // .I
                // fprintf(stderr, "    push (. vs %zu %c) +0 +1 {%zu, %zu} (%zu)\n", rhs_ref, rhs_supremal.sequence.str[rhs_idx % rhs_width], lhs_idx, rhs_idx + 1, lhs_idx * rhs_size + rhs_idx + 1);
                queue_push_back(&queue, lhs_idx * rhs_size + rhs_idx + 1, included);
            } // if
        } // if
        else if ((rhs_idx == 0 && lhs_ref < rhs_ref) || rhs_idx == rhs_size - 1)
        {
            uint8_t const lhs_value = get(lhs_dfa, lhs_idx);
            bool const lhs_deletion = lhs_value & DFA_DELETION;
            bool const lhs_insertion = lhs_value & DFA_INSERTION;
            bool const lhs_match = lhs_idx % lhs_width < lhs_width - 1 &&
                reference[lhs_supremal.start + lhs_idx / lhs_width] == lhs_supremal.sequence.str[lhs_idx % lhs_width];

            if (lhs_match)
            {
                // ..
                // fprintf(stderr, "    push (%zu MU vs .) +0 +0 {%zu, %zu} (%zu)\n", lhs_ref, lhs_idx + lhs_width + 1, rhs_idx, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx);
                queue_push_front(&queue, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx, included);
            } // if

            if (lhs_deletion)
            {
                // D.
                // fprintf(stderr, "    push (%zu DELTA vs .) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, lhs_idx + lhs_width, rhs_idx, (lhs_idx + lhs_width) * rhs_size + rhs_idx);
                queue_push_back(&queue, (lhs_idx + lhs_width) * rhs_size + rhs_idx, included);
            } // if

            if (lhs_insertion)
            {
                // I.
                // fprintf(stderr, "    push (%zu %c vs .) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, lhs_supremal.sequence.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx, (lhs_idx + 1) * rhs_size + rhs_idx);
                queue_push_back(&queue, (lhs_idx + 1) * rhs_size + rhs_idx, included);
            } // if
        } // if
        else
        {
            uint8_t const lhs_value = get(lhs_dfa, lhs_idx);
            bool const lhs_deletion = lhs_value & DFA_DELETION;
            bool const lhs_insertion = lhs_value & DFA_INSERTION;
            bool const lhs_match = lhs_idx % lhs_width < lhs_width - 1 &&
                reference[lhs_supremal.start + lhs_idx / lhs_width] == lhs_supremal.sequence.str[lhs_idx % lhs_width];

            uint8_t const rhs_value = get(rhs_dfa, rhs_idx);
            bool const rhs_deletion = rhs_value & DFA_DELETION;
            bool const rhs_insertion = rhs_value & DFA_INSERTION;
            bool const rhs_match = rhs_idx % rhs_width < rhs_width - 1 &&
                reference[rhs_supremal.start + rhs_idx / rhs_width] == rhs_supremal.sequence.str[rhs_idx % rhs_width];

            if (lhs_match && rhs_match)
            {
                // ..
                // fprintf(stderr, "    push (%zu MU) {%zu, %zu} +0 +0 (%zu)\n", lhs_ref, lhs_idx + lhs_width + 1, rhs_idx + rhs_width + 1, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx + rhs_width + 1);
                queue_push_front(&queue, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx + rhs_width + 1, included);
            } // if

            if (lhs_deletion && rhs_deletion)
            {
                // DD
                // fprintf(stderr, "    push (%zu DELTA) {%zu, %zu} +1 +0 (%zu)\n", lhs_ref, lhs_idx + lhs_width, rhs_idx + rhs_width, (lhs_idx + lhs_width) * rhs_size + rhs_idx + rhs_width);
                queue_push_front(&queue, (lhs_idx + lhs_width) * rhs_size + rhs_idx + rhs_width, included + 1);
            } // if

            if (lhs_insertion && rhs_insertion && lhs_supremal.sequence.str[lhs_idx % lhs_width] == rhs_supremal.sequence.str[rhs_idx % rhs_width])
            {
                // II
                // fprintf(stderr, "    push (%zu %c) {%zu, %zu} +1 +0 (%zu)\n", lhs_ref, lhs_supremal.sequence.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx + 1, (lhs_idx + 1) * rhs_size + rhs_idx + 1);
                queue_push_front(&queue, (lhs_idx + 1) * rhs_size + rhs_idx + 1, included + 1);
            } // if

            if (lhs_deletion && rhs_match)
            {
                // D.
                // fprintf(stderr, "    push (%zu DELTA vs %zu MU) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, rhs_ref, lhs_idx + lhs_width, rhs_idx + rhs_width + 1, (lhs_idx + lhs_width) * rhs_size + rhs_idx + rhs_width + 1);
                queue_push_back(&queue, (lhs_idx + lhs_width) * rhs_size + rhs_idx + rhs_width + 1, included);
            } // if

            if (lhs_match && rhs_deletion)
            {
                // .D
                // fprintf(stderr, "    push (%zu MU vs %zu DELTA) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, rhs_ref, lhs_idx + lhs_width + 1, rhs_idx + rhs_width, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx + rhs_width);
                queue_push_back(&queue, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx + rhs_width, included);
            } // if

            if (lhs_insertion && (rhs_deletion || rhs_match || (rhs_insertion &&
                lhs_supremal.sequence.str[lhs_idx % lhs_width] != rhs_supremal.sequence.str[rhs_idx % rhs_width])))
            {
                // I.
                // fprintf(stderr, "    push (%zu %c vs .) {%zu, %zu} +0 +1 (%zu)\n", lhs_ref, lhs_supremal.sequence.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx, (lhs_idx + 1) * rhs_size + rhs_idx);
                queue_push_back(&queue, (lhs_idx + 1) * rhs_size + rhs_idx, included);
            } // if

            if (rhs_insertion && (lhs_deletion || lhs_match || (lhs_insertion &&
                lhs_supremal.sequence.str[lhs_idx % lhs_width] != rhs_supremal.sequence.str[rhs_idx % rhs_width])))
            {
                // .I
                // fprintf(stderr, "    push (. vs %zu %c) {%zu, %zu} +0 +1 (%zu)\n", rhs_ref, rhs_supremal.sequence.str[rhs_idx % rhs_width], lhs_idx, rhs_idx + 1, lhs_idx * rhs_size + rhs_idx + 1);
                queue_push_back(&queue, lhs_idx * rhs_size + rhs_idx + 1, included);
            } // if
        } // else
    } // while
    size_t const included = queue.entries[HASH_TABLE_INDEX(queue.entries, lhs_size * rhs_size - 1)].included;

    queue_destroy(&queue);

    return included;
} // dfa_max_overlap


// FIXME: DEBUG
void
dfa_dot(GVA_Allocator const allocator, FILE* const stream,
    size_t const len, char const reference[static restrict len],
    GVA_Variant const supremal, uint8_t const dfa[static restrict 1])
{
    size_t const width = supremal.sequence.len + 1;
    size_t const size = (supremal.end - supremal.start + 1) * width;

    struct
    {
        char     seen;
        gva_uint next;
    }* queue = allocator.allocate(allocator.context, NULL, 0, sizeof(*queue) * size);
    if (queue == NULL)
    {
        return;  // OOM
    } // if

    memset(queue, 0, sizeof(*queue) * size);

    fprintf(stream, "strict digraph{\nrankdir=LR\nnode[fixedsize=true,shape=circle,width=1]\ni[label=\"\",shape=none,width=0]\ni->0\n");

    size_t head = 0;
    size_t tail = head;
    queue[head].seen = true;
    queue[head].next = GVA_NULL;
    while (head != GVA_NULL)
    {
        uint8_t const value = get(dfa, head);
        if (value & DFA_DELETION)
        {
            fprintf(stream, "%zu->%zu[label=\"%zu &Delta;\"]\n", head, head + width, supremal.start + head / width);
            if (!queue[head + width].seen)
            {
                queue[head + width].seen = true;
                queue[head + width].next = GVA_NULL;
                queue[tail].next = head + width;
                tail = head + width;
            } // if
        } // if
        if (value & DFA_INSERTION)
        {
            fprintf(stream, "%zu->%zu[label=\"%zu %c\"]\n", head, head + 1, supremal.start + head / width, supremal.sequence.str[head % width]);
            if (!queue[head + 1].seen)
            {
                queue[head + 1].seen = true;
                queue[head + 1].next = GVA_NULL;
                queue[tail].next = head + 1;
                tail = head + 1;
            } // if
        } // if
        if (head % width < width - 1 &&
            reference[supremal.start + head / width] == supremal.sequence.str[head % width])
        {
            fprintf(stream, "%zu->%zu[label=\"%zu &Mu;\"]\n", head, head + width + 1, supremal.start + head / width);
            if (!queue[head + width + 1].seen)
            {
                queue[head + width + 1].seen = true;
                queue[head + width + 1].next = GVA_NULL;
                queue[tail].next = head + width + 1;
                tail = head + width + 1;
            } // if
        } // if
        head = queue[head].next;
    } // while

    queue = allocator.allocate(allocator.context, queue, sizeof(*queue) * size, 0);

    fprintf(stream, "%zu[peripheries=2]\n}\n", size - 1);
} // dfa_dot
