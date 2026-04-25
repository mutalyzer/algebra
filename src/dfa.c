#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint8_t
#include <string.h>     // memset

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_edges
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "../include/variant.h"     // GVA_Variant
#include "array.h"              // ARRAY_*, array_*
#include "dfa.h"                // dfa*
#include "priority_queue.h"     // Priority_Queue, priority_queue_*


#include <stdio.h>      // DEBUG


static uint8_t const DFA_DELETION  = 0x1;
static uint8_t const DFA_INSERTION = 0x2;


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
                (i < end) * DFA_DELETION | (j < len) * DFA_INSERTION);
        } // for
    } // for
} // edge


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

    memset(dfas + array_length(dfas), 0, len);
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
                edge(dfas + array_length(dfas), width, variant.start + k - supremal.start, variant.end + k - supremal.start,
                    variant.sequence.len, variant.sequence.str - supremal.sequence.str + k);
            } // for
        } // for
    } // for
    array_header(dfas)->length += len;

    return dfas;
} // dfa_from_lcs_graph


size_t
dfa_max_overlap(GVA_Allocator const allocator,
    size_t const len, char const reference[static restrict len],
    GVA_Variant const lhs_supremal, uint8_t const lhs_dfa[static restrict 1],
    GVA_Variant const rhs_supremal, uint8_t const rhs_dfa[static restrict 1])
{
    size_t const lhs_width = lhs_supremal.sequence.len + 1;
    size_t const rhs_width = rhs_supremal.sequence.len + 1;
    size_t const lhs_size = (lhs_supremal.end - lhs_supremal.start + 1) * lhs_width;
    size_t const rhs_size = (rhs_supremal.end - rhs_supremal.start + 1) * rhs_width;

    Priority_Queue fringe = priority_queue_init(allocator, lhs_size * rhs_size);
    priority_queue_push(&fringe, 0, 0, 0);

    size_t count = 0;
    while (!priority_queue_empty(fringe))
    {
        count += 1;
        size_t const idx = priority_queue_pop(&fringe);

        size_t const lhs_idx = idx / rhs_size;
        size_t const rhs_idx = idx % rhs_size;

        size_t const lhs_ref = lhs_idx / lhs_width + lhs_supremal.start;
        size_t const rhs_ref = rhs_idx / rhs_width + rhs_supremal.start;

        fprintf(stderr, "{%zu, %zu} (%zu) :: %u %u\n", lhs_idx, rhs_idx, idx,
            fringe.states[idx].included, fringe.states[idx].excluded);

        if (lhs_idx == lhs_size - 1 && rhs_idx == rhs_size - 1)
        {
            fprintf(stdout, "%u %zu %zu\n", fringe.states[idx].included, count, array_length(fringe.heap));
            break;
        } // if

        if ((lhs_idx == 0 && rhs_ref < lhs_ref) || lhs_idx == lhs_size - 1)
        {
            uint8_t const value = get(rhs_dfa, rhs_idx);

            if (value & DFA_DELETION)
            {
                fprintf(stderr, "    push (. vs %zu DELTA) +0 +1 {%zu, %zu} (%zu)\n", rhs_ref, lhs_idx, rhs_idx + rhs_width, lhs_idx * rhs_size + rhs_idx + rhs_width);
                priority_queue_push(&fringe, lhs_idx * rhs_size + rhs_idx + rhs_width,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (rhs_idx % rhs_width < rhs_width - 1 &&
                reference[rhs_supremal.start + rhs_idx / rhs_width] == rhs_supremal.sequence.str[rhs_idx % rhs_width])
            {
                fprintf(stderr, "    push (. vs %zu MU) +0 +0 {%zu, %zu} (%zu)\n", rhs_ref, lhs_idx, rhs_idx + rhs_width + 1, lhs_idx * rhs_size + rhs_idx + rhs_width + 1);
                priority_queue_push(&fringe, lhs_idx * rhs_size + rhs_idx + rhs_width + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded);
            } // if

            if (value & DFA_INSERTION)
            {
                fprintf(stderr, "    push (. vs %zu %c) +0 +1 {%zu, %zu} (%zu)\n", rhs_ref, rhs_supremal.sequence.str[rhs_idx % rhs_width], lhs_idx, rhs_idx + 1, lhs_idx * rhs_size + rhs_idx + 1);
                priority_queue_push(&fringe, lhs_idx * rhs_size + rhs_idx + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if
            continue;
        } // if

        if ((rhs_idx == 0 && lhs_ref < rhs_ref) || rhs_idx == rhs_size - 1)
        {
            uint8_t const value = get(lhs_dfa, lhs_idx);

            if (value & DFA_DELETION)
            {
                fprintf(stderr, "    push (%zu DELTA vs .) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, lhs_idx + lhs_width, rhs_idx, (lhs_idx + lhs_width) * rhs_size + rhs_idx);
                priority_queue_push(&fringe, (lhs_idx + lhs_width) * rhs_size + rhs_idx,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (lhs_idx % lhs_width < lhs_width - 1 &&
                reference[lhs_supremal.start + lhs_idx / lhs_width] == lhs_supremal.sequence.str[lhs_idx % lhs_width])
            {
                fprintf(stderr, "    push (%zu MU vs .) +0 +0 {%zu, %zu} (%zu)\n", lhs_ref, lhs_idx + lhs_width + 1, rhs_idx, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx);
                priority_queue_push(&fringe, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx,
                    fringe.states[idx].included, fringe.states[idx].excluded);
            } // if

            if (value & DFA_INSERTION)
            {
                fprintf(stderr, "    push (%zu %c vs .) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, lhs_supremal.sequence.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx, (lhs_idx + 1) * rhs_size + rhs_idx);
                priority_queue_push(&fringe, (lhs_idx + 1) * rhs_size + rhs_idx,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if
            continue;
        } // if

        uint8_t const lhs_value = get(lhs_dfa, lhs_idx);
        uint8_t const rhs_value = get(rhs_dfa, rhs_idx);

        if (lhs_value & DFA_DELETION)
        {
            if (rhs_value & DFA_DELETION)
            {
                fprintf(stderr, "    push (%zu DELTA) {%zu, %zu} +1 +0 (%zu)\n", lhs_ref, lhs_idx + lhs_width, rhs_idx + rhs_width, (lhs_idx + lhs_width) * rhs_size + rhs_idx + rhs_width);
                priority_queue_push(&fringe, (lhs_idx + lhs_width) * rhs_size + rhs_idx + rhs_width,
                    fringe.states[idx].included + 1, fringe.states[idx].excluded);
            } // if

            if (rhs_idx % rhs_width < rhs_width - 1 &&
                reference[rhs_supremal.start + rhs_idx / rhs_width] == rhs_supremal.sequence.str[rhs_idx % rhs_width])
            {
                fprintf(stderr, "    push (%zu DELTA vs %zu MU) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, rhs_ref, lhs_idx + lhs_width, rhs_idx + rhs_width + 1, (lhs_idx + lhs_width) * rhs_size + rhs_idx + rhs_width + 1);
                priority_queue_push(&fringe, (lhs_idx + lhs_width) * rhs_size + rhs_idx + rhs_width + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (rhs_value & DFA_INSERTION)
            {
                fprintf(stderr, "    push (. vs %zu %c) {%zu, %zu} +0 +1 (%zu)\n", rhs_ref, rhs_supremal.sequence.str[rhs_idx % rhs_width], lhs_idx, rhs_idx + 1, lhs_idx * rhs_size + rhs_idx + 1);
                priority_queue_push(&fringe, lhs_idx * rhs_size + rhs_idx + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if
        } // if

        if (lhs_idx % lhs_width < lhs_width - 1 &&
            reference[lhs_supremal.start + lhs_idx / lhs_width] == lhs_supremal.sequence.str[lhs_idx % lhs_width])
        {
            if (rhs_value & DFA_DELETION)
            {
                fprintf(stderr, "    push (%zu MU vs %zu DELTA) +0 +1 {%zu, %zu} (%zu)\n", lhs_ref, rhs_ref, lhs_idx + lhs_width + 1, rhs_idx + rhs_width, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx + rhs_width);
                priority_queue_push(&fringe, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx + rhs_width,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (rhs_idx % rhs_width < rhs_width - 1 &&
                reference[rhs_supremal.start + rhs_idx / rhs_width] == rhs_supremal.sequence.str[rhs_idx % rhs_width])
            {
                fprintf(stderr, "    push (%zu MU) {%zu, %zu} +0 +0 (%zu)\n", lhs_ref, lhs_idx + lhs_width + 1, rhs_idx + rhs_width + 1, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx + rhs_width + 1);
                priority_queue_push(&fringe, (lhs_idx + lhs_width + 1) * rhs_size + rhs_idx + rhs_width + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded);
            } // if

            if (rhs_value & DFA_INSERTION)
            {
                fprintf(stderr, "    push (. vs %zu %c) {%zu, %zu} +0 +1 (%zu)\n", rhs_ref, rhs_supremal.sequence.str[rhs_idx % rhs_width], lhs_idx, rhs_idx + 1, lhs_idx * rhs_size + rhs_idx + 1);
                priority_queue_push(&fringe, lhs_idx * rhs_size + rhs_idx + 1,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if
        } // if

        if (lhs_value & DFA_INSERTION)
        {
            if (rhs_value & DFA_DELETION)
            {
                fprintf(stderr, "    push (%zu %c vs .) {%zu, %zu} +0 +1 (%zu)\n", lhs_ref, lhs_supremal.sequence.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx, (lhs_idx + 1) * rhs_size + rhs_idx);
                priority_queue_push(&fringe, (lhs_idx + 1) * rhs_size + rhs_idx,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (rhs_idx % rhs_width < rhs_width - 1 &&
                reference[rhs_supremal.start + rhs_idx / rhs_width] == rhs_supremal.sequence.str[rhs_idx % rhs_width])
            {
                fprintf(stderr, "    push (%zu %c vs .) {%zu, %zu} +0 +1 (%zu)\n", lhs_ref, lhs_supremal.sequence.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx, (lhs_idx + 1) * rhs_size + rhs_idx);
                priority_queue_push(&fringe, (lhs_idx + 1) * rhs_size + rhs_idx,
                    fringe.states[idx].included, fringe.states[idx].excluded + 1);
            } // if

            if (rhs_value & DFA_INSERTION)
            {
                if (lhs_supremal.sequence.str[lhs_idx % lhs_width] == rhs_supremal.sequence.str[rhs_idx % rhs_width])
                {
                    fprintf(stderr, "    push (%zu %c) {%zu, %zu} +1 +0 (%zu)\n", lhs_ref, lhs_supremal.sequence.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx + 1, (lhs_idx + 1) * rhs_size + rhs_idx + 1);
                    priority_queue_push(&fringe, (lhs_idx + 1) * rhs_size + rhs_idx + 1,
                        fringe.states[idx].included + 1, fringe.states[idx].excluded);
                } // if
                else
                {
                    fprintf(stderr, "    push (%zu %c vs .) {%zu, %zu} +0 +1 (%zu)\n", lhs_ref, lhs_supremal.sequence.str[lhs_idx % lhs_width], lhs_idx + 1, rhs_idx, (lhs_idx + 1) * rhs_size + rhs_idx);
                    priority_queue_push(&fringe, (lhs_idx + 1) * rhs_size + rhs_idx,
                        fringe.states[idx].included, fringe.states[idx].excluded + 1);

                    fprintf(stderr, "    push (. vs %zu %c) {%zu, %zu} +0 +1 (%zu)\n", rhs_ref, rhs_supremal.sequence.str[rhs_idx % rhs_width], lhs_idx, rhs_idx + 1, lhs_idx * rhs_size + rhs_idx + 1);
                    priority_queue_push(&fringe, lhs_idx * rhs_size + rhs_idx + 1,
                        fringe.states[idx].included, fringe.states[idx].excluded + 1);
                } // else
            } // if
        } // if
    } // while

    priority_queue_destroy(&fringe);

    return 0;
} // dfa_max_overlap


void
dfa_dot_traverse(size_t const len, char const reference[static len],
    GVA_Variant const supremal, uint8_t const dfa[static 1], size_t const idx)
{
    size_t const width = supremal.sequence.len + 1;
    size_t const size = (supremal.end - supremal.start + 1) * width;

    if (idx == size - 1)
    {
        return;
    } // if

    uint8_t const value = get(dfa, idx);
    if (value & DFA_DELETION)
    {
        fprintf(stderr, "%zu->%zu[label=\"%zu &Delta;\"]\n", idx, idx + width, supremal.start + idx / width);
        dfa_dot_traverse(len, reference, supremal, dfa, idx + width);
    } // if
    if (value & DFA_INSERTION)
    {
        fprintf(stderr, "%zu->%zu[label=\"%zu %c\"]\n", idx, idx + 1, supremal.start + idx / width, supremal.sequence.str[idx % width]);
        dfa_dot_traverse(len, reference, supremal, dfa, idx + 1);
    } // if
    if (reference[supremal.start + idx / width] == supremal.sequence.str[idx % width])
    {
        fprintf(stderr, "%zu->%zu[label=\"%zu &Mu;\"]\n", idx, idx + width + 1, supremal.start + idx / width);
        dfa_dot_traverse(len, reference, supremal, dfa, idx + width + 1);
    } // if
} // dfa_dot_traverse


void
dfa_dot(size_t const len, char const reference[static restrict len],
    GVA_Variant const supremal, uint8_t const dfa[static restrict 1])
{
    size_t const width = supremal.sequence.len + 1;
    size_t const size = (supremal.end - supremal.start + 1) * width;

    fprintf(stderr, "strict digraph{\nrankdir=LR\nnode[fixedsize=true,shape=circle,width=1]\ni[label=\"\",shape=none,width=0]\ni->0\n");
    dfa_dot_traverse(len, reference, supremal, dfa, 0);
    fprintf(stderr, "%zu[peripheries=2]\n}\n", size - 1);
} // dfa_dot
