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
gva_string_reverse(GVA_String self)
{
    size_t i = self.len - 1;
    size_t j = 0;
    while (i > j)
    {
        char const ch = self.str[i];
        ((char*) self.str)[i] = self.str[j];
        ((char*) self.str)[j] = ch;
        i -= 1;
        j += 1;
    } // while
} // gva_string_reverse


inline size_t
gva_string_prefix_length(GVA_String const lhs, GVA_String const rhs)
{
    size_t idx = 0;
    while (idx < lhs.len && idx < rhs.len && lhs.str[idx] == rhs.str[idx])
    {
        idx += 1;
    } // while
    return idx;
} // gva_string_prefix_length
