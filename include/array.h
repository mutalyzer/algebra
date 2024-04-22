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
    VA_Allocator const* allocator;
    size_t capacity;
    size_t item_size;
    size_t length;
} VA_Array;


void*
va_array_init(VA_Allocator const* const allocator, size_t const capacity, size_t const item_size);


void*
va_array_ensure(void* const self, size_t const count);


void*
va_array_destroy(void* const self);


#define \
va_array_append(self, value) ( \
    (self) = va_array_ensure(self, 1), \
    (self)[((VA_Array*)(self) - 1)->length] = (value), \
    &(self)[((VA_Array*)(self) - 1)->length++] \
) // va_array_append


#ifdef __cplusplus
} // extern "C"
#endif

#endif
