#ifndef ARRAY_H
#define ARRAY_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>     // NULL, size_t

#include "alloc.h"      // VA_Allocator


typedef struct
{
    size_t capacity;
    size_t length;
} VA_Array;


void*
va_array_init(VA_Allocator const allocator, size_t const capcity, size_t const item_size);


VA_Array*
va_array_header(void* const ptr);


size_t
va_array_length(void* const ptr);


void*
va_array_grow(VA_Allocator const allocator, void* const restrict ptr, size_t const capacity, size_t const item_size);


size_t
va_array_shift(void* const ptr, size_t const index, size_t const item_size);


#define va_array_destroy(allocator, ptr) ((ptr) == NULL ? NULL : allocator.alloc(allocator.context, va_array_header(ptr), va_array_header(ptr)->capacity * sizeof(*(ptr)), 0))


#define va_array_insert(allocator, ptr, value, index) do {              \
    (ptr) = va_array_grow(allocator, ptr, 1, sizeof(*ptr));             \
    if ((ptr) == NULL) break;                                           \
    (ptr)[va_array_shift(ptr, index, sizeof(*ptr))] = (value);          \
    ++va_array_header(ptr)->length;                                     \
} while (0)


#define va_array_append(allocator, ptr, value) va_array_insert(allocator, ptr, value, va_array_length(ptr))


#ifdef __cplusplus
} // extern "C"
#endif

#endif
