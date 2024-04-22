#include <stddef.h>     // NULL, size_t
#include <stdlib.h>     // EXIT_SUCCESS, free, realloc

#include "../include/alloc.h"   // VA_Allocator


#include "../src/array.c"


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


int
main(void)
{
    int* a = va_array_init(&std_allocator, 1, sizeof(*a));
    va_array_append(a, 42);
    va_array_append(a, 43);
    va_array_append(a, 44);
    a = va_array_destroy(a);

    return EXIT_SUCCESS;
} // main
