#ifndef GVA_BITSET_H
#define GVA_BITSET_H


#include <stddef.h>     // size_t


#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph


size_t*
bitset_init(GVA_Allocator const allocator, size_t const size);


void
bitset_add(size_t bitset[static 1], size_t const value);


size_t
bitset_intersection_cnt(size_t const lhs[static 1], size_t const rhs[static 1]);


size_t*
bitset_destroy(GVA_Allocator const allocator, size_t bitset[static 1]);


void
bitset_fill(GVA_LCS_Graph const graph,
            gva_uint const start, gva_uint const start_accent, gva_uint const end_accent,
            size_t dels[static restrict 1],
            size_t as[static restrict 1],
            size_t cs[static restrict 1],
            size_t gs[static restrict 1],
            size_t ts[static restrict 1]);

#endif  // GVA_BITSET_H
