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
bitset_fill(GVA_LCS_Graph graph, gva_uint start, gva_uint start_accent, gva_uint end_accent, size_t* dels, size_t* as, size_t* cs, size_t* gs, size_t* ts);

#endif  // GVA_BITSET_H
