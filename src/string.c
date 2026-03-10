#include <stddef.h>     // NULL, size_t
#include <string.h>     // memcpy

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/string.h"      // GVA_String, gva_string_*


inline GVA_String
gva_string_init(GVA_Allocator const allocator, size_t const len)
{
    GVA_String string = {len, allocator.allocate(allocator.context, NULL, 0, len)};
    if (string.str == NULL)
    {
        return (GVA_String) {0};
    } // if
    return string;
} // gva_string_init


inline void
gva_string_destroy(GVA_Allocator const allocator, GVA_String self)
{
    self.len = 0;
    self.str = allocator.allocate(allocator.context, (char*) self.str, self.len, 0);
} // gva_string_destroy


inline GVA_String
gva_string_concat(GVA_Allocator const allocator, GVA_String lhs, GVA_String const rhs)
{
    size_t const len = lhs.len + rhs.len;
    lhs.str = allocator.allocate(allocator.context, (char*) lhs.str, lhs.len, len);
    if (lhs.str == NULL)
    {
        return (GVA_String) {0};
    } // if

    memcpy((char*) lhs.str + lhs.len, rhs.str, rhs.len);
    lhs.len = len;
    return lhs;
} // gva_string_concat


inline GVA_String
gva_string_dup(GVA_Allocator const allocator, GVA_String const self)
{
    return gva_string_concat(allocator, (GVA_String) {0}, self);
} // gva_string_dup


inline void
gva_string_reverse(size_t const len, char str[static len])
{
    size_t i = 0;
    size_t j = len - 1;
    while (i < j)
    {
        char const ch = str[i];
        str[j] = str[i];
        str[i] = ch;
        i += 1;
        j -= 1;
    } // while
} // gva_string_reverse


inline size_t
gva_prefix_length(size_t const len_lhs, char const lhs[static restrict len_lhs],
    size_t const len_rhs, char const rhs[static restrict len_rhs])
{
    size_t i = 0;
    while (i < len_lhs && i < len_rhs && lhs[i] == rhs[i])
    {
        i += 1;
    } // while
    return i;
} // gva_prefix_length
