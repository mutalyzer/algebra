#ifndef BITSET_H
#define BITSET_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdbool.h>    // bool
#include <stddef.h>     // size_t

#include "alloc.h"      // VA_Allocator


typedef struct VA_Bitset VA_Bitset;


VA_Bitset*
va_bitset_init(VA_Allocator const allocator, size_t const capacity);


VA_Bitset*
va_bitset_destroy(VA_Allocator const allocator, VA_Bitset* const bitset);


VA_Bitset*
va_bitset_set(VA_Bitset* const bitset, size_t const index);


bool
va_bitset_test(VA_Bitset const* const bitset, size_t const index);


#ifdef __cplusplus
} // extern "C"
#endif

#endif
