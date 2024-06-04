#include <stddef.h>     // NULL, size_t

#include "../include/alloc.h"   // VA_Allocator
#include "../include/array.h"   // VA_Array, va_array_*


void*
va_array_init(VA_Allocator const allocator[static 1], size_t const capacity, size_t const item_size)
{
    size_t const min_capacity = capacity < 4 ? 4 : capacity;
    VA_Array* const header = allocator->alloc(allocator->context, NULL, 0, sizeof(*header) + min_capacity * item_size);  // OVERFLOW
    if (header == NULL)
    {
        return NULL;
    } // if
    header->capacity = min_capacity;
    header->item_size = item_size;
    header->length = 0;
    return header + 1;
} // va_array_init


void*
va_array_ensure(VA_Allocator const allocator[static restrict 1], void* const restrict ptr, size_t const count)
{
    VA_Array* header = (VA_Array*)(ptr) - 1;
    if (header->capacity - header->length >= count)
    {
        return ptr;
    } // if

    size_t const old_size = sizeof(*header) + header->capacity * header->item_size;
    header->capacity *= 2;  // OVERFLOW
    while (header->capacity - header->length < count)
    {
        header->capacity *= 2;  // OVERFLOW
    } // while

    header = allocator->alloc(allocator->context, header, old_size, sizeof(*header) + header->capacity * header->item_size);  // OVERFLOW
    if (header == NULL)
    {
        return NULL;
    } // if
    return header + 1;
} // va_array_ensure


void*
va_array_destroy(VA_Allocator const allocator[static restrict 1], void* const restrict ptr)
{
    if (ptr == NULL)
    {
        return NULL;
    } // if
    VA_Array* const header = (VA_Array*)(ptr) - 1;
    return allocator->alloc(allocator->context, header, sizeof(*header) + header->capacity * header->item_size, 0);
} // va_array_destroy
