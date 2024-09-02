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
va_std_alloc(void* const restrict context, void* const restrict ptr, size_t const old_size, size_t const size)
{
    (void) context;
    (void) old_size;

    if (size == 0)
    {
        free(ptr);
        return NULL;
    } // if

    char* const restrict new_ptr = realloc(ptr, size);
    if (new_ptr == NULL)
    {
        free(ptr);
        return NULL;  // OOM
    } // if

    return new_ptr;
} // va_std_alloc


static VA_Allocator const va_std_allocator = { .alloc = va_std_alloc };


#ifdef __cplusplus
} // extern "C"
#endif

#endif
