#ifndef STATIC_ALLOC_H
#define STATIC_ALLOC_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>     // NULL, size_t


#define KiB(size) ((size_t) size * 1024ULL)
#define MiB(size) (KiB(size) * 1024ULL)
#define GiB(size) (MiB(size) * 1024ULL)


#define VA_STATIC_BUFFER_SIZE MiB(50)


typedef struct
{
    size_t _size;
    size_t free;
    char data[VA_STATIC_BUFFER_SIZE];
} VA_Static_Buffer;


static VA_Static_Buffer va_static_buffer = { ._size = sizeof(va_static_buffer) };


static void*
va_static_alloc(void* const restrict context, void* const restrict ptr, size_t const old_size, size_t const size)
{
    VA_Static_Buffer* const restrict buffer = context;

    if (size == 0)
    {
        return NULL;
    } // if

    if (old_size >= size)
    {
        return ptr;
    } // if

    if (ptr != NULL && old_size > 0)
    {
        if (sizeof(buffer->data) - buffer->free < old_size)
        {
            return NULL;  // OOM
        } // if
        char const* const restrict src = ptr;
        for (size_t i = 0; i < old_size; ++i)
        {
            buffer->data[buffer->free + i] = src[i];
        } // for
    } // if

    if (sizeof(buffer->data) - buffer->free < size)
    {
        return NULL;  // OOM
    } // if

    void* const restrict new_ptr = buffer->data + buffer->free;
    buffer->free += size;
    return new_ptr;
} // va_static_alloc


static VA_Allocator const va_static_allocator =
{
    va_static_alloc,
    &va_static_buffer
};


#ifdef __cplusplus
} // extern "C"
#endif

#endif
