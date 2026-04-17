#include <stddef.h>     // size_t
#include <stdio.h>      // stderr, fprintf  FIXME: DEBUG

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_edges
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "../include/variant.h"     // GVA_Variant
#include "array.h"      // ARRAY_*, array_length
#include "nfa.h"        // NFA*, nfa_*


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
