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
    VA_Bitset* const self = allocator.alloc(allocator.context, NULL, 0, sizeof(*self) + slots * size);  // OVERFLOW
    if (self == NULL)
    {
        return NULL;  // OOM
    } // if
    self->capacity = capacity;
    for (size_t i = 0; i < slots; ++i)
    {
        self->data[i] = 0;
    } // for
    return self;
} // va_bitset_init


inline VA_Bitset*
va_bitset_destroy(VA_Allocator const allocator, VA_Bitset* const self)
{
    if (self == NULL)
    {
        return NULL;
    } // if
    size_t const size = sizeof(self->data[0]) * CHAR_BIT;
    size_t const slots = (self->capacity + size - 1) / size;
    return allocator.alloc(allocator.context, self, sizeof(*self) + slots * size, 0);
} // va_bitset_destroy


inline void
va_bitset_set(VA_Bitset* const self, size_t const index)
{
    if (self == NULL || self->capacity <= index)
    {
        return;
    } // if
    size_t const size = sizeof(self->data[0]) * CHAR_BIT;
    self->data[index / size] |= ((size_t) 1) << (index % size);
} // va_bitset_set


inline bool
va_bitset_test(VA_Bitset const* const self, size_t const index)
{
    if (self == NULL || self->capacity <= index)
    {
        return false;
    } // if
    size_t const size = sizeof(self->data[0]) * CHAR_BIT;
    return self->data[index / size] & (((size_t) 1) << (index % size));
} // va_bitset_test
