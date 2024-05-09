#ifndef STD_ALLOC_H
#define STD_ALLOC_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>     // NULL, size_t
#include <stdlib.h>     // free, realloc
#include <string.h>     // memset

#include "alloc.h"      // VA_Allocator


static inline void*
va_std_alloc(void* const context, void* const ptr, size_t const old_size, size_t const size)
{
    (void) context;
    (void) old_size;
    if (size == 0)
    {
        free(ptr);
        return NULL;
    } // if
    char* const here = realloc(ptr, size);
    if (size > old_size)
    {
        (void) memset(here + old_size, 0, size - old_size);
    } // if
    return here;
} // va_std_alloc


static VA_Allocator const va_std_allocator = { .alloc = va_std_alloc };


#ifdef __cplusplus
} // extern "C"
#endif

#endif
