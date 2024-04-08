#include <stdbool.h>    // bool
#include <stddef.h>     // size_t

#include "../include/variant.h"     // VA_*, va_*


bool
va_variant_eq(VA_Variant const lhs, VA_Variant const rhs)
{
    return lhs.start == rhs.start && lhs.end == rhs.end &&
           lhs.obs_start == rhs.obs_start && lhs.obs_end == rhs.obs_end;
} // va_variant_eq


size_t
va_variant_len(VA_Variant const variant)
{
    return variant.end - variant.start + variant.obs_end - variant.obs_start;
} // va_variant_len
