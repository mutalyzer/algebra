#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // printf
#include <stdlib.h>     // EXIT_*

#include "../include/alloc.h"       // VA_Allocator
#include "../include/array.h"       // VA_Array


#define KiB(size) ((size_t) size * 1024ULL)
#define MiB(size) (KiB(size) * 1024ULL)
#define GiB(size) (MiB(size) * 1024ULL)


typedef struct
{
    size_t _size;
    size_t free;
    char buffer[MiB(50)];
} Static_Buffer;


static Static_Buffer static_buffer = { ._size = sizeof(static_buffer) };


static void*
static_alloc(void* const context, void* const ptr, size_t const old_size, size_t const size)
{
    Static_Buffer* const buffer = (Static_Buffer*) context;

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
        if (sizeof(buffer->buffer) - buffer->free < old_size)
        {
            return NULL;  // OOM
        } // if
        char const* const src = ptr;
        for (size_t i = 0; i < old_size; ++i)
        {
            buffer->buffer[buffer->free + i] = src[i];
        } // for
    } // if

    if (sizeof(buffer->buffer) - buffer->free < size)
    {
        return NULL;  // OOM
    } // if

    void* const here = buffer->buffer + buffer->free;
    buffer->free += size;
    return here;
} // static_alloc


static VA_Allocator const static_allocator =
{
    static_alloc,
    &static_buffer
};


int
main(int argc, char* argv[argc + 1])
{
    (void) argv;


    int* a = va_array_init(&static_allocator, 100, sizeof(*a));

    a[0] = 42;

    printf("%p\n", (void*) a);
    printf("%d\n", a[0]);
    printf("%d\n", a[1]);

    a = va_array_destroy(&static_allocator, a);

    return EXIT_SUCCESS;
} // main
