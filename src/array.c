#include <stddef.h>     // NULL, size_t


#include "../include/allocator.h"   // GVA_Allocator
#include "array.h"       // GVA_Array, gva_array_*


inline void*
array_init(GVA_Allocator const allocator, size_t const capacity,
    size_t const item_size)
{
    if (capacity == 0)
    {
        return NULL;  // Empty array
    } // if
    Array* const header = allocator.allocate(
        allocator.context,
        NULL,
        0,
        sizeof(*header) + capacity * item_size  // OVERFLOW
    );
    if (header == NULL)
    {
        return NULL;  // OOM
    } // if
    header->capacity = capacity;
    header->length = 0;
    return header + 1;
} // array_init


inline size_t
array_length(void const* const self)
{
    return self == NULL ? 0 : array_header(self)->length;
} // array_length


inline Array*
array_header(void const* const self)
{
    return (Array*) self - 1;
} // array_header


void*
array_ensure(GVA_Allocator const allocator, void* const self,
    size_t const item_size, size_t const extra)
{
    // Initial allocation
    if (self == NULL)
    {
        return array_init(allocator, extra, item_size);
    } // if

    // The current capacity is sufficient.
    Array* const restrict header = array_header(self);
    if (header->capacity - header->length >= extra)
    {
        return header + 1;
    } // if

    // A simple geometric growth reallocation scheme.
    size_t new_capacity = header->capacity * 2;  // OVERFLOW
    while (new_capacity - header->length < extra)
    {
        new_capacity *= 2;  // OVERFLOW
    } // while

    Array* const restrict new_header = allocator.allocate(
        allocator.context,
        header,
        sizeof(*header) + header->capacity * item_size,
        sizeof(*new_header) + new_capacity * item_size  // OVERFLOW
    );
    if (new_header == NULL)
    {
        return NULL;  // OOM
    } // if

    new_header->capacity = new_capacity;
    return new_header + 1;
} // array_ensure_one
