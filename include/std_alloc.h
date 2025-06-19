// NOT FREESTANDING
#ifndef GVA_STD_ALLOC_H
#define GVA_STD_ALLOC_H


#include <stddef.h>         // NULL, size_t
#include <stdlib.h>         // free, realloc

#include "allocator.h"      // GVA_Allocator


// Provides a simple global allocation strategy based on `realloc` from
// libc.
static inline void*
gva_std_allocate(void* const restrict context, void const* const restrict ptr, size_t const old_size, size_t const new_size)
{
    (void) context;
    (void) old_size;

    if (new_size == 0)
    {
        free((void*) ptr);
        return NULL;
    } // if

    void* const new_ptr = realloc((void*) ptr, new_size);
    if (new_ptr == NULL)
    {
        free((void*) ptr);
        return NULL;  // OOM
    } // if

    return new_ptr;
} // va_std_allocate


static GVA_Allocator const gva_std_allocator = { .allocate = gva_std_allocate };


#endif // GVA_STD_ALLOC_H
