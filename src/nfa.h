#ifndef GVA_NFA_H
#define GVA_NFA_H


#include <stddef.h>     // size_t
#include <stdint.h>     // uint8_t

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph
#include "../include/string.h"      // GVA_String
#include "../include/types.h"       // gva_uint


typedef struct
{
    uint8_t deletion  : 1;
    uint8_t insertion : 1;
    uint8_t match     : 1;
} DFA_State;  // FIXME: could realy be 3 bits
// Thm. each state has a maximum of two outgoing transitions.


typedef struct
{
    size_t     size;
    GVA_String observed;
    DFA_State* states;
} DFA;


DFA
dfa_from_lcs_graph(GVA_Allocator const allocator, GVA_LCS_Graph const graph);


void
dfa_destroy(GVA_Allocator const allocator, DFA self[static 1]);


size_t
dfa_max_overlap(GVA_Allocator const allocator, DFA const lhs, DFA const rhs);


// FIXME: DEBUG
void
dfa_dot(DFA const self);


static char const NFA_LAMBDA =   0;
static char const NFA_DELETE = 127;


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
