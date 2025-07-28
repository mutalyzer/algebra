#ifndef GVA_BITSET_H
#define GVA_BITSET_H


#include <stddef.h>     // size_t

#include "../include/allocator.h"   // GVA_Allocator


size_t*
bitset_init(GVA_Allocator const allocator, size_t const size);


void
bitset_add(size_t self[static 1], size_t const start, size_t const end);


size_t
bitset_intersection_cnt(size_t const lhs[static 1], size_t const rhs[static 1]);


size_t*
bitset_destroy(GVA_Allocator const allocator, size_t const* const self);


#endif  // GVA_BITSET_H
