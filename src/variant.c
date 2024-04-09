#include <stdbool.h>    // bool
#include <stddef.h>     // size_t

#include "../include/variant.h"     // VA_*, va_*


static void
string_copy(size_t const len, char const source[static len],
            char dest[static len])
{
    for (size_t i = 0; i < len; ++i)
    {
        dest[i] = source[i];
    } // for
} // string_copy


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


size_t
va_patch(size_t const len_ref, char const reference[static len_ref],
         size_t const len_obs, char const observed[static len_obs],
         size_t const len_var, VA_Variant const variants[len_var],
         size_t const len, char string[static len])
{
    size_t idx = 0;
    size_t start = 0;
    for (size_t i = 0; i < len_var; ++i)
    {
        size_t len_slice = variants[i].start - start;
        if (len - idx < len_slice)
        {
            return 0;  // not enough space
        } // if
        string_copy(len_slice, reference + start, string + idx);
        idx += len_slice;

        len_slice = variants[i].obs_end - variants[i].obs_start;
        if (len - idx < len_slice)
        {
            return 0;  // not enough space
        } // if
        string_copy(len_slice, observed + variants[i].obs_start, string + idx);
        idx += len_slice;

        start = variants[i].end;
    } // for

    size_t const len_slice = len_ref - start;
    if (len - idx < len_slice)
    {
        return 0;  // not enough space
    } // if
    string_copy(len_slice, reference + start, string + idx);
    return idx + len_slice;
} // va_patch
