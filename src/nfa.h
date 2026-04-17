#ifndef GVA_NFA_H
#define GVA_NFA_H


#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph
#include "../include/types.h"       // gva_uint


static char const NFA_DELETE = 127;
static char const NFA_LAMBDA = (char) 235;


typedef struct
{
    char     symbol;
    gva_uint tail;
    gva_uint next;
} NFA_Transition;


typedef struct
{
    gva_uint*       states;
    NFA_Transition* transitions;
    gva_uint        source;
    gva_uint        sink;
} NFA;


NFA
nfa_from_lcs_graph(GVA_Allocator const allocator, GVA_LCS_Graph const graph);


void
nfa_destroy(GVA_Allocator const allocator, NFA self[static 1]);


// FIXME: DEBUG
void
nfa_dot(NFA const self);


#endif  // GVA_NFA_H
