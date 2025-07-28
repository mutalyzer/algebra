#ifndef GVA_ALLOCATOR_H
#define GVA_ALLOCATOR_H


#include <stddef.h>     // size_t


// Generic allocator interface. Provides a single `allocate` function
// similar to `realloc` from libc. The following behavior is expected:
// - when `new_size` is zero, the function behaves like `free` and
//   returns `NULL`.
// - when `new_size` is not zero, the function behaves like `realloc`.
//   Returns `NULL` if the allocation request cannot be fulfilled.
//
// An opaque `context` may be provided and allows for custom allocators
// that can be scoped/threaded/aligned/etc.
//
// See `std_alloc.h` for a simple libc wrapper.
typedef struct
{
    void* (*allocate)(void* const restrict context, void* const restrict ptr, size_t const old_size, size_t const new_size);
    void* context;
} GVA_Allocator;


#endif // GVA_ALLOCATOR_H
