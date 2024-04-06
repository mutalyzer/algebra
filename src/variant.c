#include <stddef.h>     // size_t

#include "../include/variant.h"     // VA_*, va_*


size_t
va_variant_len(VA_Variant const variant)
{
    return variant.end - variant.start + variant.obs_end - variant.obs_start;
} // va_variant_len
