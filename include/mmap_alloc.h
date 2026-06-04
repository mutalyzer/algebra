#ifndef GVA_MMAP_ALLOC_H
#define GVA_MMAP_ALLOC_H


#include <stddef.h>     // size_t


typedef struct
{
    void*  addr;  // this *must* be the first field
    size_t len;
    int    fd;
} GVA_Mmap_Context;


// Provides a simple file backed memory mapped allocation strategy.
void*
gva_mmap_allocate(void* const restrict context, void* const restrict ptr,
    size_t const old_size, size_t const new_size);


GVA_Mmap_Context
gva_mmap_context_init(char const path[static 1], size_t const initial_len);


void
gva_mmap_context_destroy(GVA_Mmap_Context const ctx[static 1]);


#endif // GVA_MMAP_ALLOC_H
