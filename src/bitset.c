// FIXME: __builtin_popcountll
#include <limits.h>     // CHAR_BIT
#include <stddef.h>     // NULL, size_t


#include "../include/allocator.h"   // GVA_Allocator
#include "array.h"  // ARRAY_*, array_*


inline size_t*
bitset_init(GVA_Allocator const allocator, size_t const size)
{
    size_t const len = (size + sizeof(size_t) * CHAR_BIT - 1) / (sizeof(size_t) * CHAR_BIT);
    size_t* bitset = array_init(allocator, len, sizeof(*bitset));
    if (bitset == NULL)
    {
        return NULL;
    } // if
    array_header(bitset)->length = len;
    for (size_t i = 0; i < len; ++i)
    {
        bitset[i] = 0;
    } // for
    return bitset;
} // bitset_init


inline void
bitset_add(size_t bitset[static 1], size_t const value)
{
    bitset[value / (sizeof(size_t) * CHAR_BIT)] |= (1ULL << (value % (sizeof(size_t) * CHAR_BIT)));
} // bitset_add


inline size_t
bitset_intersection_cnt(size_t const lhs[static 1], size_t const rhs[static 1])
{
    size_t const len = array_length(lhs);
    size_t count = 0;
    for (size_t i = 0; i < len; ++i)
    {
        count += __builtin_popcountll(lhs[i] & rhs[i]);
    } // for
    return count;
} // bitset_intersection_cnt


inline size_t*
bitset_destroy(GVA_Allocator const allocator, size_t bitset[static 1])
{
    return ARRAY_DESTROY(allocator, bitset);
} // bitset_destroy
