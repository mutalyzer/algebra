#ifndef GVA_DFA_H
#define GVA_DFA_H


#include <stdbool.h>    // bool
#include <stddef.h>     // size_t
#include <stdint.h>     // uint8_t

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph
#include "../include/variant.h"     // GVA_Variant


#include <stdio.h>      // FIXME: DEBUG


uint8_t*
dfa_init(GVA_Allocator const allocator,
    size_t const height, size_t const width);


uint8_t*
dfa_concat(GVA_Variant const lhs, uint8_t lhs_dfa[static restrict 1],
    GVA_Variant const rhs, uint8_t const rhs_dfa[static restrict 1],
    size_t const offset);


uint8_t*
dfa_from_alignment(GVA_Allocator const allocator, uint8_t* restrict dfas,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs]);


uint8_t*
dfa_from_lcs_graph(GVA_Allocator const allocator, uint8_t* dfas,
    GVA_LCS_Graph const graph, size_t const start, size_t const end);


size_t
dfa_overlap(GVA_Allocator const allocator,
    size_t const len, char const reference[static restrict len],
    GVA_Variant const lhs, uint8_t const lhs_dfa[static restrict 1],
    GVA_Variant const rhs, uint8_t const rhs_dfa[static restrict 1]);


bool
dfa_disjoint(GVA_Variant const lhs, uint8_t const lhs_dfa[static restrict 1],
    GVA_Variant const rhs, uint8_t const rhs_dfa[static restrict 1]);


// FIXME: DEBUG
void
dfa_dot(FILE* restrict const stream,
    size_t const len, char const reference[static restrict len],
    GVA_Variant const variant, uint8_t const dfa[static restrict 1]);


void
dfa_svg(FILE* restrict const stream,
    size_t const len, char const reference[static restrict len],
    GVA_Variant const variant, uint8_t const dfa[static restrict 1]);


#endif  // GVA_DFA_H
