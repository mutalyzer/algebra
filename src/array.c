#include <stddef.h>     // NULL, size_t

#include "../include/alloc.h"   // VA_Allocator
#include "../include/array.h"   // VA_Array, va_array_*


inline void*
va_array_init(VA_Allocator const allocator, size_t const capacity, size_t const item_size)
{
    return va_array_grow(allocator, NULL, capacity, item_size);
} // va_array_init


inline VA_Array*
va_array_header(void* const ptr)
{
    return (VA_Array*) (ptr) - 1;
} // va_array_header


inline size_t
va_array_length(void* const ptr)
{
    return ptr == NULL ? 0 : va_array_header(ptr)->length;
} // va_array_length


void*
va_array_grow(VA_Allocator const allocator, void* const restrict ptr, size_t const capacity, size_t const item_size)
{
    if (capacity == 0)
    {
        return ptr;
    } // if

    if (ptr == NULL)
    {
        VA_Array* const restrict header = allocator.alloc(allocator.context, NULL, 0, sizeof(*header) + capacity * item_size);  // OVERFLOW
        if (header == NULL)
        {
            return NULL;  // OOM
        } // if
        header->capacity = capacity;
        header->length = 0;
        return header + 1;
    } // if

    VA_Array* restrict header = va_array_header(ptr);
    if (header->capacity - header->length >= capacity)
    {
        return ptr;
    } // if

    size_t new_capacity = header->capacity * 2;  // OVERFLOW
    while (new_capacity < capacity)
    {
        new_capacity *= 2;  // OVERFLOW
    } // while

    header = allocator.alloc(allocator.context, header, sizeof(*header) + header->capacity * item_size, sizeof(*header) * new_capacity * item_size);
    if (header == NULL)
    {
        return NULL;  // OOM
    } // if

    header->capacity = new_capacity;
    return header + 1;
} // va_array_grow


size_t
va_array_shift(void* const ptr, size_t const index, size_t const item_size)
{
    size_t const length = va_array_length(ptr);
    if (index >= length)
    {
        return length;
    } // if

    for (size_t i = 0; i < (length - index) * item_size; ++i)
    {
        ((char*) ptr)[length * item_size + item_size - 1 - i] = ((char*) ptr)[length * item_size - 1 - i];
    } // for
    return index;
} // va_array_shift
