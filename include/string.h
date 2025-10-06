#ifndef GVA_STRING_H
#define GVA_STRING_H


#include <stddef.h>     // size_t

#include "allocator.h"  // GVA_Allocator


#define GVA_STRING_PRINT(string) (int) string.len, string.str
#define GVA_STRING_FMT "%.*s"


typedef struct
{
    size_t      len;
    char const* str;
} GVA_String;


void
gva_string_destroy(GVA_Allocator const allocator, GVA_String self);


GVA_String
gva_string_concat(GVA_Allocator const allocator,
    GVA_String lhs, GVA_String const rhs);


GVA_String
gva_string_dup(GVA_Allocator const allocator, GVA_String const self);


#endif
