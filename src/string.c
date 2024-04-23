#include <stddef.h>     // NULL, size_t

#include "../include/alloc.h"       // VA_Allocator
#include "../include/string.h"      // VA_String, va_string_*


VA_String
va_string_concat(VA_Allocator const* const allocator, VA_String lhs, VA_String const rhs)
{
    lhs.data = allocator->alloc(allocator->context, lhs.data, lhs.length, lhs.length + rhs.length);  // OVERFLOW
    if (lhs.data == NULL)
    {
        return (VA_String) {0, NULL};
    } // if

    for (size_t i = 0; i < rhs.length; ++i)
    {
        lhs.data[lhs.length + i] = rhs.data[i];
    } // for
    lhs.length += rhs.length;
    return lhs;
} // va_string_concat


VA_String
va_string_slice(VA_Allocator const* const allocator, VA_String const string, size_t const start, size_t const end)
{
    (void) allocator;
    (void) string;
    (void) start;
    (void) end;
    return (VA_String) {0, NULL};
} // va_string_slice
