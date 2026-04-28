#ifndef GVA_DFA_H
#define GVA_DFA_H


#include <stddef.h>     // size_t
#include <stdint.h>     // uint8_t

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph
#include "../include/variant.h"     // GVA_Variant


uint8_t*
dfa_from_alignment(GVA_Allocator const allocator, uint8_t* dfas,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs]);


uint8_t*
dfa_from_lcs_graph(GVA_Allocator const allocator, uint8_t* dfas,
    GVA_LCS_Graph const graph);


size_t
dfa_max_overlap(GVA_Allocator const allocator,
    size_t const len, char const reference[static restrict len],
    GVA_Variant const lhs_supremal, uint8_t const lhs_dfa[static restrict 1],
    GVA_Variant const rhs_supremal, uint8_t const rhs_dfa[static restrict 1]);


// FIXME: DEBUG
void
dfa_dot(size_t const len, char const reference[static restrict len],
    GVA_Variant const supremal, uint8_t const dfa[static restrict 1]);


#endif  // GVA_DFA_H
