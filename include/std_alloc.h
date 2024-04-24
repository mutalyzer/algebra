#ifndef STD_ALLOC_H
#define STD_ALLOC_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>     // NULL, size_t
#include <stdlib.h>     // free, realloc

#include "alloc.h"      // VA_Allocator


static inline void*
std_alloc(void* const context, void* const ptr, size_t const old_size, size_t const size)
{
    (void) context;
    (void) old_size;
    if (size == 0)
    {
        free(ptr);
        return NULL;
    } // if
    return realloc(ptr, size);
} // std_alloc


static VA_Allocator const std_allocator = { .alloc = std_alloc };


#ifdef __cplusplus
} // extern "C"
#endif

#endif
