#include <stdbool.h>    // bool
#include <stddef.h>     // NULL, size_t
#include <string.h>     // memcmp, memcpy

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/string.h"      // GVA_String, gva_string_*, gva_prefix_length
#include "../include/variant.h"     // gva_patch, gva_variant_*,
#include "common.h"     // MIN


inline bool
gva_variant_eq(GVA_Variant const lhs, GVA_Variant const rhs)
{
    return lhs.start == rhs.start && lhs.end == rhs.end &&
        lhs.sequence.len == rhs.sequence.len &&
        memcmp(lhs.sequence.str, rhs.sequence.str, MIN(lhs.sequence.len, rhs.sequence.len)) == 0;
} // gva_variant_eq


inline GVA_Variant
gva_variant_dup(GVA_Allocator const allocator, GVA_Variant const variant)
{
    return (GVA_Variant) {variant.start, variant.end, gva_string_dup(allocator, variant.sequence)};
} // gva_variant_dup


inline size_t
gva_variant_length(GVA_Variant const variant)
{
    return variant.end - variant.start + variant.sequence.len;
} // gva_variant_length


inline GVA_Variant
gva_variant_prefix_trimmed(size_t const len_ref, char const reference[static len_ref],
    GVA_Variant const variant)
{
    size_t const len = gva_prefix_length(variant.end - variant.start, reference + variant.start,
        variant.sequence.len, variant.sequence.str);
    return (GVA_Variant) {variant.start + len, variant.end, {variant.sequence.len - len, variant.sequence.str + len}};
} // gva_variant_prefix_trimmed


GVA_String
gva_patch(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const n, GVA_Variant const variants[static restrict n])
{
    size_t deleted = 0;
    size_t inserted = 0;
    for (size_t i = 0; i < n; ++i)
    {
        deleted += variants[i].end - variants[i].start;
        inserted += variants[i].sequence.len;
    } // for

    char* str = allocator.allocate(allocator.context, NULL, 0, len_ref + inserted - deleted);
    if (str == NULL)
    {
        return (GVA_String) {0};
    } // if

    size_t len = 0;
    size_t start = 0;
    for (size_t i = 0; i < n; ++i)
    {
        memcpy(str + len, reference + start, variants[i].start - start);
        len += variants[i].start - start;
        memcpy(str + len, variants[i].sequence.str, variants[i].sequence.len);
        len += variants[i].sequence.len;
        start = variants[i].end;
    } // for

    if (start < len_ref)
    {
        memcpy(str + len, reference + start, len_ref - start);
        len += len_ref - start;
    } // if

    return (GVA_String) {len, str};
} // gva_patch
