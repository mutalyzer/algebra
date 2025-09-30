#include <stddef.h>     // NULL, size_t
#include <string.h>     // memcpy

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/string.h"      // GVA_String, gva_string_*


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
        return (GVA_String) {0, NULL};
    } // if

    memcpy((char*) lhs.str + lhs.len, rhs.str, rhs.len);
    lhs.len = len;
    return lhs;
} // gva_string_concat


inline GVA_String
gva_string_dup(GVA_Allocator const allocator, GVA_String const self)
{
    return gva_string_concat(allocator, (GVA_String) {0, NULL}, self);
} // gva_string_dup
