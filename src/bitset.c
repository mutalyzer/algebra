#include <limits.h>     // CHAR_BIT
#include <stdbool.h>    // bool, false
#include <stddef.h>     // NULL, size_t

#include "../include/alloc.h"   // VA_Allocator
#include "../include/bitset.h"  // va_bitset_*


struct VA_Bitset
{
    size_t capacity;
    size_t data[];
};


inline VA_Bitset*
va_bitset_init(VA_Allocator const allocator, size_t const capacity)
{
    size_t const size = sizeof((VA_Bitset) {0}.data[0]) * CHAR_BIT;
    size_t const slots = (capacity + size - 1) / size;
    VA_Bitset* const bitset = allocator.alloc(allocator.context, NULL, 0, sizeof(*bitset) + slots * size);  // OVERFLOW
    if (bitset == NULL)
    {
        return NULL;  // OOM
    } // if
    bitset->capacity = capacity;
    for (size_t i = 0; i < slots; ++i)
    {
        bitset->data[i] = 0;
    } // for
    return bitset;
} // va_bitset_init


inline VA_Bitset*
va_bitset_destroy(VA_Allocator const allocator, VA_Bitset* const bitset)
{
    if (bitset == NULL)
    {
        return NULL;
    } // if
    size_t const size = sizeof(bitset->data[0]) * CHAR_BIT;
    size_t const slots = (bitset->capacity + size - 1) / size;
    return allocator.alloc(allocator.context, bitset, sizeof(*bitset) + slots * size, 0);
} // va_bitset_destroy


inline VA_Bitset*
va_bitset_set(VA_Bitset* const bitset, size_t const index)
{
    if (bitset == NULL || bitset->capacity <= index)
    {
        return bitset;
    } // if
    size_t const size = sizeof(bitset->data[0]) * CHAR_BIT;
    bitset->data[index / size] |= ((size_t) 1) << (index % size);
    return bitset;
} // va_bitset_set


inline bool
va_bitset_test(VA_Bitset const* const bitset, size_t const index)
{
    if (bitset == NULL || bitset->capacity <= index)
    {
        return false;
    } // if
    size_t const size = sizeof(bitset->data[0]) * CHAR_BIT;
    return bitset->data[index / size] & (((size_t) 1) << (index % size));
} // va_bitset_test
