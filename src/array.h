#ifndef GVA_ARRAY_H
#define GVA_ARRAY_H


#include <stddef.h>         // NULL, size_t

#include "../include/allocator.h"      // GVA_Allocator


// Internal: the metadata of a dynamic array.
// Should not be used directly.
typedef struct
{
    size_t capacity;
    size_t length;
} Array;


// Creates a new dynamic array allocated with `capacity` number of
// initial elements.
// It is the responsibility of the caller to provide the same allocator
// to any consequent allocation events (including the destruction of the
// array).
// Usage:
//     int* a = array_init(gva_std_allocator, 10, sizeof(*a));
// When `capacity` is zero the result is `NULL`. Any `NULL` value behaves
// as an empty array.
void*
array_init(GVA_Allocator const allocator, size_t const capacity,
    size_t const item_size);


// Returns the number of elements currently in use by the array.
size_t
array_length(void* const self);


// Internal: Returns the metadata of the array.
// `self` must not be `NULL`.
Array*
array_header(void* const self);


// Internal: Returns the possibly (re)allocated array with place for
// one additional element.
void*
array_ensure_one(GVA_Allocator const allocator, void* const self,
    size_t const item_size);


// Destroys an array. Does not destroy the elements of the array.
// Returns `NULL`.
// Usage:
//     a = ARRAY_DESTROY(gva_std_allocator, a);
// Warning: `self` is evaluated multiple times.
#define ARRAY_DESTROY(allocator, self) (                                \
    (self) == NULL ? NULL :                                             \
        allocator.allocate(                                             \
            allocator.context,                                          \
            array_header(self),                                         \
            array_header(self)->capacity * sizeof(*(self)),             \
            0                                                           \
        )                                                               \
)


// Adds an element at the end of the array. Returns the new length of the
// array. 0 if (re)allocation failed.
// Warning: `self` is evaluated multiple times.
#define ARRAY_APPEND(allocator, self, value) (                          \
    (self) = array_ensure_one(allocator, (self), sizeof(*(self))),      \
    (self) == NULL ? 0 : (                                              \
        (self)[array_header(self)->length] = (value),                   \
        array_header(self)->length += 1                                 \
    )                                                                   \
)


#endif // GVA_ARRAY_H
