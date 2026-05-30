#ifndef GVA_STD_ALLOC_H
#define GVA_STD_ALLOC_H


#include <stddef.h>         // size_t

#include "allocator.h"      // GVA_Allocator


// Provides a simple global allocation strategy based on `realloc` from libc.
void*
gva_std_allocate(void* const restrict context, void* const restrict ptr,
    size_t const old_size, size_t const new_size);


// Global general purpose allocator.
static GVA_Allocator const gva_std_allocator = { .allocate = gva_std_allocate };


#endif // GVA_STD_ALLOC_H
