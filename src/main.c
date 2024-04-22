#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // stderr, printf
#include <stdlib.h>     // EXIT_*, free, realloc

#include "../include/alloc.h"   // VA_Allocator
#include "../include/array.h"   // va_array_*


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
main(int argc, char* argv[argc + 1])
{
    (void) argv;

    int* a = va_array_init(&std_allocator, 1, sizeof(*a));

    va_array_append(a, 42);
    va_array_append(a, 43);

    printf("%d\n", a[1]);

    a = va_array_destroy(a);

    return EXIT_SUCCESS;
} // main
