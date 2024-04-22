#include <stddef.h>     // NULL, size_t

#include "../include/alloc.h"   // VA_Allocator
#include "../include/array.h"   // VA_Array, va_array_*


void*
va_array_init(VA_Allocator const* const allocator, size_t const capacity, size_t const item_size)
{
    VA_Array* const header = allocator->alloc(allocator->context, NULL, 0, sizeof(*header) + capacity * item_size);  // OVERFLOW
    if (header == NULL)
    {
        return NULL;
    } // if
    header->allocator = allocator;
    header->capacity = capacity;
    header->item_size = item_size;
    header->length = 0;
    return header + 1;
} // va_array_init


void*
va_array_ensure(void* const self, size_t const count)
{
    VA_Array* header = (VA_Array*)(self) - 1;
    size_t const free = header->capacity - header->length;
    if (free >= count)
    {
        return self;
    } // if

    size_t const capacity = count - free > 2 * header->capacity ? count - free + 16 / header->item_size : 2 * header->capacity;  // OVERFLOW
    header = header->allocator->alloc(header->allocator->context, header, sizeof(*header) + header->capacity * header->item_size, sizeof(*header) + capacity * header->item_size);
    if (header == NULL)
    {
        return NULL;
    } // if
    header->capacity += capacity;
    return header + 1;
} // va_array_ensure


void*
va_array_destroy(void* const self)
{
    VA_Array* const header = (VA_Array*)(self) - 1;
    return header->allocator->alloc(header->allocator->context, header, sizeof(*header) + header->capacity * header->item_size, 0);
} // va_array_destroy
