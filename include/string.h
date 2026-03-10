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


GVA_String
gva_string_init(GVA_Allocator const allocator, size_t const len);


void
gva_string_destroy(GVA_Allocator const allocator, GVA_String self);


GVA_String
gva_string_concat(GVA_Allocator const allocator,
    GVA_String lhs, GVA_String const rhs);


GVA_String
gva_string_dup(GVA_Allocator const allocator, GVA_String const self);


void
gva_string_reverse(size_t const len, char str[static len]);


size_t
gva_prefix_length(size_t const len_lhs, char const lhs[static restrict len_lhs],
    size_t const len_rhs, char const rhs[static restrict len_rhs]);


#endif // GVA_STRING_H
