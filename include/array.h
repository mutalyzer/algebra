#ifndef ARRAY_H
#define ARRAY_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>     // size_t

#include "alloc.h"      // VA_Allocator


typedef struct
{
    size_t capacity;
    size_t item_size;
    size_t length;
} VA_Array;


void*
va_array_init(VA_Allocator const* const allocator, size_t const capacity, size_t const item_size);


void*
va_array_ensure(VA_Allocator const* const allocator, void* const ptr, size_t const count);


void*
va_array_destroy(VA_Allocator const* const allocator, void* const ptr);


#define \
va_array_append(allocator, ptr, value) (                    \
    (ptr) = va_array_ensure(allocator, ptr, 1),             \
    (ptr)[((VA_Array*)(ptr) - 1)->length] = (value),        \
    &(ptr)[((VA_Array*)(ptr) - 1)->length++]                \
) // va_array_append


#ifdef __cplusplus
} // extern "C"
#endif

#endif
