#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // stderr, fprintf  FIXME: DEBUG
#include <string.h>     // memset

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_edges
#include "../include/string.h"      // GVA_String, gva_string_*
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "../include/variant.h"     // GVA_Variant
#include "array.h"              // ARRAY_*, array_length
#include "nfa.h"                // DFA, dfa*, NFA*, nfa_*
#include "priority_queue.h"     // Priority_Queue, priority_queue_*


static inline void
set(uint8_t dfa[static 1], size_t const idx, uint8_t const value)
{
    dfa[idx / 4] |= value << (2 * (idx % 4));
} // set


static inline uint8_t
get(uint8_t const dfa[static 1], size_t const idx)
{
    return (dfa[idx / 4] >> (2 * (idx % 4)) & 0x3);
} // get


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
                (i < end) * DFA2_DELETION | (j < len) * DFA2_INSERTION);
        } // for
    } // for
} // edge


uint8_t*
dfa2_from_lcs_graph(GVA_Allocator const allocator,
    uint8_t* dfa, GVA_LCS_Graph const graph)
{
    GVA_Variant const supremal = gva_lcs_graph_supremal(graph);
    size_t const width = supremal.sequence.len + 1;
    size_t const size = (supremal.end - supremal.start + 1) * width;
    dfa = array_ensure(allocator, dfa, 1, (size + 1) / 4);
    if (dfa == NULL)
    {
        return NULL;  // OOM
    } // if

    array_header(dfa)->length += (size + 1) / 4;
    memset(dfa, 0, array_length(dfa));

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
                edge(dfa, width, variant.start + k - supremal.start, variant.end + k - supremal.start,
                    variant.sequence.len, variant.sequence.str - graph.observed.str + k);
            } // for
        } // for
    } // for

    return dfa;
} // dfa2_from_lcs_graph


void
dfa2_dot_traverse(size_t const len, char const reference[static len],
    GVA_Variant const supremal, uint8_t const dfa[static 1], size_t const idx)
{
    size_t const width = supremal.sequence.len + 1;
    size_t const size = (supremal.end - supremal.start + 1) * width;

    if (idx == size - 1)
    {
        return;
    } // if

    uint8_t const value = get(dfa, idx);
    if (value & DFA2_DELETION)
    {
        fprintf(stderr, "%zu->%zu[label=\"%zu &Delta;\"]\n", idx, idx + width, supremal.start + idx / width);
        dfa2_dot_traverse(len, reference, supremal, dfa, idx + width);
    } // if
    if (value & DFA2_INSERTION)
    {
        fprintf(stderr, "%zu->%zu[label=\"%zu %c\"]\n", idx, idx + 1, supremal.start + idx / width, supremal.sequence.str[idx % width]);
        dfa2_dot_traverse(len, reference, supremal, dfa, idx + 1);
    } // if
    if (reference[idx / width] == supremal.sequence.str[idx % width])
    {
        fprintf(stderr, "%zu->%zu[label=\"%zu &Mu;\"]\n", idx, idx + width + 1, supremal.start + idx / width);
        dfa2_dot_traverse(len, reference, supremal, dfa, idx + width + 1);
    } // if
} // dfa2_dot_traverse


void
dfa2_dot(size_t const len, char const reference[static len],
    GVA_Variant const supremal, uint8_t const dfa[static 1])
{
    size_t const width = supremal.sequence.len + 1;
    size_t const size = (supremal.end - supremal.start + 1) * width;

    fprintf(stderr, "strict digraph{\nrankdir=LR\nnode[fixedsize=true,shape=circle,width=1]\ni[label=\"\",shape=none,width=0]\ni->0\n");
    dfa2_dot_traverse(len, reference, supremal, dfa, 0);
    fprintf(stderr, "%zu[peripheries=2]\n}\n", size - 1);
} // dfa2_dot


static inline void
dfa_edge(DFA self[static 1], size_t const width,
    size_t const start, size_t const end,
    size_t const len, size_t const offset)
{
    for (size_t i = start; i <= end; ++i)
    {
        for (size_t j = 0; j <= len; ++j)
        {
            if (i < end)
            {
                self->states[i * width + offset + j].deletion = 1;
            } // if
            if (j < len)
            {
                self->states[i * width + offset + j].insertion = 1;
            } // if
        } // for
    } // for
} // dfa_edge


DFA
dfa_from_lcs_graph(GVA_Allocator const allocator, GVA_LCS_Graph const graph)
{
    GVA_Variant const supremal = gva_lcs_graph_supremal(graph);
    size_t const size = (supremal.end - supremal.start + 1) * (supremal.sequence.len + 1);
    DFA dfa =
    {
        .start = supremal.start,
        .observed = supremal.sequence,
        .states = allocator.allocate(allocator.context, NULL, 0, size),
    };
    if (dfa.states == NULL)
    {
        dfa.size = 0;
        return dfa;  // OOM
    } // if

    dfa.size = size;
    memset(dfa.states, 0, dfa.size);

    size_t const width = dfa.observed.len + 1;
    for (size_t i = 0; i < array_length(graph.nodes); ++i)
    {
        for (size_t j = 0; j < graph.nodes[i].match.length; ++j)
        {
            dfa.states[(graph.nodes[i].match.row + j - dfa.start) * width + graph.nodes[i].match.col + j - dfa.start].match = 1;
        } // for
        for (gva_uint j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
        {
            GVA_Variant variant = {0};
            size_t const count = gva_edges(graph.observed.str,
                graph.nodes[i].match, graph.nodes[graph.edges[j].tail].match,
                i == graph.source, graph.nodes[graph.edges[j].tail].edges == GVA_NULL,
                &variant);
            for (size_t k = 0; k < count; ++k)
            {
                dfa_edge(&dfa, width, variant.start + k - dfa.start, variant.end + k - dfa.start,
                    variant.sequence.len, variant.sequence.str - dfa.observed.str + k);
            } // for
        } // for
    } // for
    return dfa;
} // dfa_from_lcs_graph


inline void
dfa_destroy(GVA_Allocator const allocator, DFA self[static 1])
{
    self->states = allocator.allocate(allocator.context, self->states, self->size * sizeof(*self->states), 0);
    self->size = 0;
} // dfa_destroy


size_t
dfa_max_overlap(GVA_Allocator const allocator, DFA const lhs, DFA const rhs)
{
    size_t const lhs_width = lhs.observed.len + 1;
    size_t const rhs_width = rhs.observed.len + 1;

    Priority_Queue fringe = priority_queue_init(allocator, lhs.size * rhs.size);
    priority_queue_push(&fringe, 0, 0, 0);

    size_t count = 0;
    while (!priority_queue_empty(fringe))
    {
        count += 1;
        size_t const idx = priority_queue_pop(&fringe);

        size_t const lhs_idx = idx / rhs.size;
        size_t const rhs_idx = idx % rhs.size;

        size_t const lhs_ref = lhs_idx / lhs_width + lhs.start;
        size_t const rhs_ref = rhs_idx / rhs_width + rhs.start;

        fprintf(stderr, "{%zu, %zu} (%zu) :: %u %u\n", lhs_idx, rhs_idx, idx,
            fringe.states[idx].included, fringe.states[idx].excluded);

        if (lhs_idx == lhs.size - 1 && rhs_idx == rhs.size - 1)
        {
            fprintf(stdout, "%u %zu %zu\n", fringe.states[idx].included, count, array_length(fringe.heap));
            break;
        } // if

        if (lhs_idx == 0 && rhs_ref < lhs_ref)
        {
            if (rhs.states[rhs_idx].deletion)
            {
                fprintf(stderr, "    push (. vs %zu DELTA) +0 +1 {%zu, %zu} (%zu)\n", rhs_ref, lhs_idx, rhs_idx + rhs_width, lhs_idx * rhs.size + rhs_idx + rhs_width);
                priority_queue_push(&fringe, lhs_idx * rhs.size + rhs_idx + rhs_width,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (rhs.states[rhs_idx].match)
            {
                fprintf(stderr, "    push (. vs %zu MU) +0 +0 {%zu, %zu} (%zu)\n", rhs_ref, lhs_idx, rhs_idx + rhs_width + 1, lhs_idx * rhs.size + rhs_idx + rhs_width + 1);
                priority_queue_push(&fringe, lhs_idx * rhs.size + rhs_idx + rhs_width + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded);
            } // if

            if (rhs.states[rhs_idx].insertion)
            {
                fprintf(stderr, "    push (. vs %zu %c) +0 +1 {%zu, %zu} (%zu)\n", rhs_ref, rhs.observed.str[rhs_idx % rhs_width], lhs_idx, rhs_idx + 1, lhs_idx * rhs.size + rhs_idx + 1);
                priority_queue_push(&fringe, lhs_idx * rhs.size + rhs_idx + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if
            continue;
        } // if

        if (rhs_idx == 0 && lhs_ref < rhs_ref)
        {
            if (lhs.states[lhs_idx].deletion)
            {
                fprintf(stderr, "    push (%zu DELTA vs .) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, lhs_idx + lhs_width, rhs_idx, (lhs_idx + lhs_width) * rhs.size + rhs_idx);
                priority_queue_push(&fringe, (lhs_idx + lhs_width) * rhs.size + rhs_idx,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (lhs.states[lhs_idx].match)
            {
                fprintf(stderr, "    push (%zu MU vs .) +0 +0 {%zu, %zu} (%zu)\n", lhs_ref, lhs_idx + lhs_width + 1, rhs_idx, (lhs_idx + lhs_width + 1) * rhs.size + rhs_idx);
                priority_queue_push(&fringe, (lhs_idx + lhs_width + 1) * rhs.size + rhs_idx,
                    fringe.states[idx].included, fringe.states[idx].excluded);
            } // if

            if (lhs.states[lhs_idx].insertion)
            {
                fprintf(stderr, "    push (%zu %c vs .) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, lhs.observed.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx, (lhs_idx + 1) * rhs.size + rhs_idx);
                priority_queue_push(&fringe, (lhs_idx + 1) * rhs.size + rhs_idx,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if
            continue;
        } // if

        if (lhs.states[lhs_idx].deletion)
        {
            if (rhs.states[rhs_idx].deletion)
            {
                fprintf(stderr, "    push (%zu DELTA) {%zu, %zu} +1 +0 (%zu)\n", lhs_ref, lhs_idx + lhs_width, rhs_idx + rhs_width, (lhs_idx + lhs_width) * rhs.size + rhs_idx + rhs_width);
                priority_queue_push(&fringe, (lhs_idx + lhs_width) * rhs.size + rhs_idx + rhs_width,
                    fringe.states[idx].included + 1, fringe.states[idx].excluded);
            } // if

            if (rhs.states[rhs_idx].match)
            {
                fprintf(stderr, "    push (%zu DELTA vs %zu MU) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, rhs_ref, lhs_idx + lhs_width, rhs_idx + rhs_width + 1, (lhs_idx + lhs_width) * rhs.size + rhs_idx + rhs_width + 1);
                priority_queue_push(&fringe, (lhs_idx + lhs_width) * rhs.size + rhs_idx + rhs_width + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (rhs.states[rhs_idx].insertion)
            {
                fprintf(stderr, "    push (. vs %zu %c) {%zu, %zu} +0 +1 (%zu)\n", rhs_ref, rhs.observed.str[rhs_idx % rhs_width], lhs_idx, rhs_idx + 1, lhs_idx * rhs.size + rhs_idx + 1);
                priority_queue_push(&fringe, lhs_idx * rhs.size + rhs_idx + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (rhs_idx == rhs.size - 1)
            {
                fprintf(stderr, "    push (%zu DELTA vs .) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, lhs_idx + lhs_width, rhs_idx, (lhs_idx + lhs_width) * rhs.size + rhs_idx);
                priority_queue_push(&fringe, (lhs_idx + lhs_width) * rhs.size + rhs_idx,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if
        } // if

        if (lhs.states[lhs_idx].match)
        {
            if (rhs.states[rhs_idx].deletion)
            {
                fprintf(stderr, "    push (%zu MU vs %zu DELTA) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, rhs_ref, lhs_idx + lhs_width + 1, rhs_idx + rhs_width, (lhs_idx + lhs_width + 1) * rhs.size + rhs_idx + rhs_width);
                priority_queue_push(&fringe, (lhs_idx + lhs_width + 1) * rhs.size + rhs_idx + rhs_width,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (rhs.states[rhs_idx].match)
            {
                fprintf(stderr, "    push (%zu MU) {%zu, %zu} +0 +0 (%zu)\n", lhs_ref, lhs_idx + lhs_width + 1, rhs_idx + rhs_width + 1, (lhs_idx + lhs_width + 1) * rhs.size + rhs_idx + rhs_width + 1);
                priority_queue_push(&fringe, (lhs_idx + lhs_width + 1) * rhs.size + rhs_idx + rhs_width + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded);
            } // if

            if (rhs.states[rhs_idx].insertion)
            {
                fprintf(stderr, "    push (. vs %zu %c) {%zu, %zu} +0 +1 (%zu)\n", rhs_ref, rhs.observed.str[rhs_idx % rhs_width], lhs_idx, rhs_idx + 1, lhs_idx * rhs.size + rhs_idx + 1);
                priority_queue_push(&fringe, lhs_idx * rhs.size + rhs_idx + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (rhs_idx == rhs.size - 1)
            {
                fprintf(stderr, "    push (%zu MU vs .) +0 +0 {%zu, %zu} (%zu)\n", lhs_ref, lhs_idx + lhs_width + 1, rhs_idx, (lhs_idx + lhs_width + 1) * rhs.size + rhs_idx);
                priority_queue_push(&fringe, (lhs_idx + lhs_width + 1) * rhs.size + rhs_idx,
                    fringe.states[idx].included, fringe.states[idx].excluded);
            } // if
        } // if

        if (lhs.states[lhs_idx].insertion)
        {
            if (rhs.states[rhs_idx].deletion)
            {
                fprintf(stderr, "    push (%zu %c vs .) {%zu, %zu} +0 +1 (%zu)\n", lhs_ref, lhs.observed.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx, (lhs_idx + 1) * rhs.size + rhs_idx);
                priority_queue_push(&fringe, (lhs_idx + 1) * rhs.size + rhs_idx,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (rhs.states[rhs_idx].match)
            {
                fprintf(stderr, "    push (%zu %c vs .) {%zu, %zu} +0 +1 (%zu)\n", lhs_ref, lhs.observed.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx, (lhs_idx + 1) * rhs.size + rhs_idx);
                priority_queue_push(&fringe, (lhs_idx + 1) * rhs.size + rhs_idx,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (rhs.states[rhs_idx].insertion)
            {
                if (lhs.observed.str[lhs_idx % lhs_width] == rhs.observed.str[rhs_idx % rhs_width])
                {
                    fprintf(stderr, "    push (%zu %c) {%zu, %zu} +1 +0 (%zu)\n", lhs_ref, lhs.observed.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx + 1, (lhs_idx + 1) * rhs.size + rhs_idx + 1);
                    priority_queue_push(&fringe, (lhs_idx + 1) * rhs.size + rhs_idx + 1,
                        fringe.states[idx].included + 1, fringe.states[idx].excluded);
                } // if
                else
                {
                    fprintf(stderr, "    push (%zu %c vs .) {%zu, %zu} +0 +1 (%zu)\n", lhs_ref, lhs.observed.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx, (lhs_idx + 1) * rhs.size + rhs_idx);
                    priority_queue_push(&fringe, (lhs_idx + 1) * rhs.size + rhs_idx,
                        fringe.states[idx].included, fringe.states[idx].excluded + 1);

                    fprintf(stderr, "    push (. vs %zu %c) {%zu, %zu} +0 +1 (%zu)\n", rhs_ref, rhs.observed.str[rhs_idx % rhs_width], lhs_idx, rhs_idx + 1, lhs_idx * rhs.size + rhs_idx + 1);
                    priority_queue_push(&fringe, lhs_idx * rhs.size + rhs_idx + 1,
                        fringe.states[idx].included, fringe.states[idx].excluded + 1);
                } // else
            } // if

            if (rhs_idx == rhs.size - 1)
            {
                fprintf(stderr, "    push (%zu %c vs .) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, lhs.observed.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx, (lhs_idx + 1) * rhs.size + rhs_idx);
                priority_queue_push(&fringe, (lhs_idx + 1) * rhs.size + rhs_idx,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if
        } // if

        if (lhs_idx == lhs.size - 1)
        {
            if (rhs.states[rhs_idx].deletion)
            {
                fprintf(stderr, "    push (. vs %zu DELTA) +0 +1 {%zu, %zu} (%zu)\n", rhs_ref, lhs_idx, rhs_idx + rhs_width, lhs_idx * rhs.size + rhs_idx + rhs_width);
                priority_queue_push(&fringe, lhs_idx * rhs.size + rhs_idx + rhs_width,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (rhs.states[rhs_idx].match)
            {
                fprintf(stderr, "    push (. vs %zu MU) +0 +0 {%zu, %zu} (%zu)\n", rhs_ref, lhs_idx, rhs_idx + rhs_width + 1, lhs_idx * rhs.size + rhs_idx + rhs_width + 1);
                priority_queue_push(&fringe, lhs_idx * rhs.size + rhs_idx + rhs_width + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded);
            } // if

            if (rhs.states[rhs_idx].insertion)
            {
                fprintf(stderr, "    push (. vs %zu %c) +0 +1 {%zu, %zu} (%zu)\n", rhs_ref, rhs.observed.str[rhs_idx % rhs_width], lhs_idx, rhs_idx + 1, lhs_idx * rhs.size + rhs_idx + 1);
                priority_queue_push(&fringe, lhs_idx * rhs.size + rhs_idx + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if
        } // if
    } // while

    priority_queue_destroy(&fringe);

    return 0;
} // dfa_max_overlap


void
dfa_dot(DFA const self)
{
    size_t const width = self.observed.len + 1;
    fprintf(stderr, "strict digraph{\nrankdir=LR\nnode[fixedsize=true,shape=circle,width=1]\ni[label=\"\",shape=none,width=0]\ni->0\n");
    for (size_t i = 0; i < self.size; ++i)
    {
        if (self.states[i].deletion)
        {
            fprintf(stderr, "%zu->%zu[label=\"%zu &Delta;\"]\n", i, i + width, self.start + i / width);
        } // if
        if (self.states[i].insertion)
        {
            fprintf(stderr, "%zu->%zu[label=\"%zu %c\"]\n", i, i + 1, self.start + i / width, self.observed.str[i % width]);
        } // if
        if (self.states[i].match)
        {
            fprintf(stderr, "%zu->%zu[label=\"%zu &Mu;\"]\n", i, i + width + 1, self.start + i / width);
        } // if
    } // for
    fprintf(stderr, "%zu[peripheries=2]\n}\n", self.size - 1);
} // dfa_dot


static void
nfa_edge(GVA_Allocator const allocator,
    NFA self[static 1], GVA_Variant const variant)
{
    for (size_t i = variant.start; i <= variant.end; ++i)
    {
        for (size_t j = 0; j <= variant.sequence.len; ++j)
        {
            size_t const tail = ARRAY_APPEND(allocator, self->states, GVA_NULL);
            if (i > variant.start)
            {
                size_t const head = tail - variant.sequence.len - 1;
                self->states[head] = ARRAY_APPEND(allocator, self->transitions,
                    ((NFA_Transition)
                    {
                        .symbol = NFA_DELETE,
                        .tail = tail,
                        .next = self->states[head],
                    }));
            } // if
            if (j > 0)
            {
                size_t const head = tail - 1;
                self->states[head] = ARRAY_APPEND(allocator, self->transitions,
                    ((NFA_Transition)
                    {
                        .symbol = variant.sequence.str[j - 1],
                        .tail = tail,
                        .next = self->states[head],
                    }));
            } // if
        } // for
    } // for
} // nfa_edge


NFA
nfa_from_lcs_graph(GVA_Allocator const allocator, GVA_LCS_Graph const graph)
{
    NFA nfa = {.source = graph.source};
    for (size_t i = 0; i < array_length(graph.nodes); ++i)
    {
        size_t const idx = ARRAY_APPEND(allocator, nfa.states, GVA_NULL);
        if (graph.nodes[i].edges == GVA_NULL)
        {
            nfa.sink = idx;
        } // if
    } // for

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
                nfa.states[i] = ARRAY_APPEND(allocator, nfa.transitions,
                    ((NFA_Transition)
                    {
                        .symbol = NFA_LAMBDA,
                        .tail = array_length(nfa.states),
                        .next = nfa.states[i],
                    }));
                nfa_edge(allocator, &nfa,
                    (GVA_Variant)
                    {
                        .start = variant.start + k,
                        .end = variant.end + k,
                        .sequence.len = variant.sequence.len,
                        .sequence.str = variant.sequence.str + k,
                    });
                size_t const head = array_length(nfa.states) - 1;
                nfa.states[head] = ARRAY_APPEND(allocator, nfa.transitions,
                    ((NFA_Transition)
                    {
                        .symbol = NFA_LAMBDA,
                        .tail = graph.edges[j].tail,
                        .next = nfa.states[head],
                    }));
            } // for
        } // for
    } // for

    return nfa;
} // nfa_from_lcs_graph


inline void
nfa_destroy(GVA_Allocator const allocator, NFA self[static 1])
{
    self->transitions = ARRAY_DESTROY(allocator, self->transitions);
    self->states = ARRAY_DESTROY(allocator, self->states);
} // nfa_destroy


void
nfa_dot(NFA const self)
{
    fprintf(stderr, "strict digraph{\nrankdir=LR\nnode[fixedsize=true,label=\"\",shape=circle,width=1]\ni[shape=none,width=0]\ni->%u\n", self.source);
    for (size_t i = 0; i < array_length(self.states); ++i)
    {
        for (gva_uint j = self.states[i]; j != GVA_NULL; j = self.transitions[j].next)
        {
            if (self.transitions[j].symbol == NFA_DELETE)
            {
                fprintf(stderr, "%zu->%u[label=\"&Delta;\"]\n", i, self.transitions[j].tail);
            } // if
            else if (self.transitions[j].symbol == NFA_LAMBDA)
            {
                fprintf(stderr, "%zu->%u[label=\"&lambda;\"]\n", i, self.transitions[j].tail);
            } // if
            else
            {
                fprintf(stderr, "%zu->%u[label=\"%c\"]\n", i, self.transitions[j].tail, self.transitions[j].symbol);
            } // else
        } // for
    } // for
    fprintf(stderr, "%u[peripheries=2]\n}\n", self.sink);
} // nfa_dot
