#ifndef STRING_H
#define STRING_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>     // size_t

#include "alloc.h"      // VA_Allocator


typedef struct
{
    size_t length;
    char* data;
} VA_String;


VA_String
va_string_concat(VA_Allocator const allocator[static 1], VA_String lhs, VA_String const rhs);


VA_String
va_string_slice(VA_Allocator const allocator[static 1], VA_String const string, size_t const start, size_t const end);


#ifdef __cplusplus
} // extern "C"
#endif

#endif
