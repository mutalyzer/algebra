#include <stdbool.h>    // bool
#include <stddef.h>     // NULL, size_t
#include <string.h>     // memcmp

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/types.h"       // gva_uint
#include "../include/string.h"      // GVA_String, gva_string_*
#include "../include/variant.h"     // gva_parse_spdi, gva_patch
#include "common.h"     // MIN


static inline bool
is_digit(char const ch)
{
    return ch >= '0' && ch <= '9';
} // is_digit


static inline bool
is_dna_nucleotide(char const ch)
{
    return ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T';
} // is_dna_nucleotide


static size_t
match_number(size_t const len, char const expression[static restrict len],
    gva_uint* const restrict num)
{
    size_t idx = 0;
    while (idx < len && is_digit(expression[idx]))
    {
        *num = *num * 10 + expression[idx] - '0';  // OVERFLOW
        idx += 1;
    } // while
    return idx;
} // match_number


static size_t
match_sequence(size_t const len, char const expression[static len])
{
    size_t idx = 0;
    while (idx < len && is_dna_nucleotide(expression[idx]))
    {
        idx += 1;
    } // while
    return idx;
} // match_sequence


static size_t
match_until(size_t const len, char const expression[static len], char const delim)
{
    size_t idx = 0;
    while (idx < len && expression[idx] != delim)
    {
        idx += 1;
    } // while
    return idx < len && expression[idx] == delim ? idx : 0;
} // match_until


size_t
gva_parse_spdi(size_t const len, char const expression[static restrict len],
    GVA_Variant variants[static restrict 1])
{
    variants[0] = (GVA_Variant) {0, 0, {0, NULL}};
    // sequence
    size_t idx = match_until(len, expression, ':');
    if (idx >= len)
    {
        return 0;  // expected ':'
    } // if
    idx += 1;  // ':'

    size_t tok = match_number(len - idx, expression + idx, &variants[0].start);
    if (tok == 0)
    {
        return 0;  // expected number (location)
    } // if
    idx += tok;

    if (idx >= len || expression[idx] != ':')
    {
        return 0;  // expected ':'
    } // if
    idx += 1;

    tok = match_number(len - idx, expression + idx, &variants[0].end);
    if (tok == 0)
    {
        tok = match_sequence(len - idx, expression + idx);
        variants[0].end = tok;  // OVERFLOW
    } // if
    idx += tok;
    variants[0].end += variants[0].start;

    if (idx >= len || expression[idx] != ':')
    {
        return 0;  // expected ':'
    } // if
    idx += 1;

    variants[0].sequence.str = expression + idx;
    tok = match_sequence(len - idx, expression + idx);
    variants[0].sequence.len = tok;
    idx += tok;

    return idx == len;
} // gva_parse_spdi


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


GVA_String
gva_patch(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const n, GVA_Variant const variants[static restrict n])
{
    GVA_String observed = {0, NULL};
    size_t start = 0;
    for (size_t i = 0; i < n; ++i)
    {
        observed = gva_string_concat(allocator, observed, (GVA_String) {variants[i].start - start, reference + start});
        observed = gva_string_concat(allocator, observed, variants[i].sequence);
        start = variants[i].end;
    } // for

    if (start < len_ref)
    {
        observed = gva_string_concat(allocator, observed, (GVA_String) {len_ref - start, reference + start});
    } // if

    return observed;
} // gva_patch
