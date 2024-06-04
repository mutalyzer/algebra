#ifndef ALLOC_H
#define ALLOC_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>     // size_t


typedef struct
{
    void* (*alloc)(void* const restrict context, void* const restrict ptr, size_t const old_size, size_t const size);
    void* context;
} VA_Allocator;


#ifdef __cplusplus
} // extern "C"
#endif

#endif
