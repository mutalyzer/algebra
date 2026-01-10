// FIXME: __builtin_popcountll
#include <limits.h>     // CHAR_BIT
#include <stddef.h>     // NULL, size_t
#include <string.h>     // memset

#include "../include/allocator.h"   // GVA_Allocator
#include "array.h"      // ARRAY_*, array_*
#include "bitset.h"     // bitset_*
#include "common.h"     // MIN


inline size_t*
bitset_init(GVA_Allocator const allocator, size_t const size)
{
    size_t const len = (size + sizeof(size_t) * CHAR_BIT - 1) / (sizeof(size_t) * CHAR_BIT);
    size_t* const bitset = array_init(allocator, len, sizeof(*bitset));
    if (bitset == NULL)
    {
        return NULL;
    } // if
    array_header(bitset)->length = len;
    return memset(bitset, 0, len * sizeof(*bitset));
} // bitset_init


inline void
bitset_add(size_t self[static 1], size_t const start, size_t const end)
{
    size_t const start_idx = start / (sizeof(*self) * CHAR_BIT);
    size_t const end_idx = (end - 1) / (sizeof(*self) * CHAR_BIT);
    for (size_t i = start_idx; i <= end_idx; ++i)
    {
        size_t mask = ~0x0ULL;
        if (i == start_idx)
        {
            mask &= ~0x0ULL << (start % (sizeof(*self) * CHAR_BIT));
        } // if
        if (i == end_idx)
        {
            mask &= ~0x0ULL >> (sizeof(*self) * CHAR_BIT - (end % (sizeof(*self) * CHAR_BIT)));
        } // if
        self[i] |= mask;
    } // for
} // bitset_add


inline size_t
bitset_intersection_cnt(size_t const lhs[static 1], size_t const rhs[static 1])
{
    size_t count = 0;
    size_t const len = MIN(array_length(lhs), array_length(rhs));
    for (size_t i = 0; i < len; ++i)
    {
        count += __builtin_popcountll(lhs[i] & rhs[i]);
    } // for
    return count;
} // bitset_intersection_cnt


inline size_t*
bitset_destroy(GVA_Allocator const allocator, size_t const* const self)
{
    return ARRAY_DESTROY(allocator, self);
} // bitset_destroy
