#include <stdbool.h>    // bool
#include <stddef.h>     // size_t

#include "../include/variant.h"     // VA_*, va_*


static void
string_copy(size_t const len, char const source[static restrict len],
            char dest[static restrict len])
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
va_patch(size_t const len_ref, char const reference[static restrict len_ref],
         size_t const len_obs, char const observed[static restrict len_obs],
         size_t const len_var, VA_Variant const variants[len_var],
         size_t const len, char string[static restrict len])
{
    size_t idx = 0;
    size_t start = 0;
    for (size_t i = 0; i < len_var; ++i)
    {
        size_t len_slice = variants[i].start - start;
        if (len - idx < len_slice)
        {
            return 0;  // not enough space in dest
        } // if
        if (len_ref - start < len_slice)
        {
            return 0;  // not in reference
        } // if
        string_copy(len_slice, reference + start, string + idx);
        idx += len_slice;

        len_slice = variants[i].obs_end - variants[i].obs_start;
        printf("%zu\n", len_slice);
        printf("%zu\n", idx);
        if (len - idx < len_slice)
        {
            return 0;  // not enough space in dest
        } // if
        if (len_obs - variants[i].obs_start < len_slice)
        {
            return 0;  // not in observed
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
