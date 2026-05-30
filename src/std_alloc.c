// NOT FREESTANDING
#include <stddef.h>     // NULL, size_t
#include <stdlib.h>     // free, realloc

#include "../include/std_alloc.h"   // gva_std_allocate


inline void*
gva_std_allocate(void* const restrict context, void* const restrict ptr,
    size_t const old_size, size_t const new_size)
{
    (void) context;   // ignored
    (void) old_size;  // ignored

    if (new_size == 0)
    {
        free(ptr);
        return NULL;
    } // if

    void* const new_ptr = realloc(ptr, new_size);
    if (new_ptr == NULL)
    {
        free(ptr);
        return NULL;  // OOM
    } // if

    return new_ptr;
} // gva_std_allocate
