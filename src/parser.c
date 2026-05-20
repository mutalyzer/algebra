#include <stdbool.h>    // bool
#include <stddef.h>     // size_t
#include <string.h>     // memcmp

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/parser.h"      // gva_parse_*
#include "../include/string.h"      // GVA_String
#include "../include/types.h"       // gva_uint
#include "../include/variant.h"     // GVA_Variant


#include <stdio.h>      // FIXME: DEBUG


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


static inline size_t
match(size_t const len, char const expression[static restrict len],
    char const word[static restrict 1])
{
    size_t idx = 0;
    while (idx < len && expression[idx] == word[idx])
    {
        idx += 1;
    } // while
    return word[idx] == '\0' ? idx : 0;
} // match


static inline size_t
match_number(size_t const len, char const expression[static restrict len],
    gva_uint value[static restrict 1])
{
    size_t idx = 0;
    while (idx < len && is_digit(expression[idx]))
    {
        *value = *value * 10 + expression[idx] - '0';  // OVERFLOW
        idx += 1;
    } // while
    return idx;
} // match_number


static inline size_t
match_location(size_t const len, char const expression[static restrict len],
    gva_uint start[static restrict 1], gva_uint end[static restrict 1])
{
    size_t idx = match_number(len, expression, start);
    if (idx == 0)
    {
        return 0;  // expected number
    } // if

    *end = *start;
    *start -= 1;
    if (idx < len && expression[idx] == '_')
    {
        idx += 1;  // '_'
        *end = 0;
        size_t const tok = match_number(len - idx, expression + idx, end);
        if (tok == 0)
        {
            return 0;  // expected number
        } // if
        idx += tok;
    } // if

    return idx;
} // match_location


static inline size_t
match_sequence(size_t const len, char const expression[static len])
{
    size_t idx = 0;
    while (idx < len && is_dna_nucleotide(expression[idx]))
    {
        idx += 1;
    } // while
    return idx;
} // match_sequence


static inline size_t
match_until(size_t const len, char const expression[static len], char const delim)
{
    size_t idx = 0;
    while (idx < len && expression[idx] != delim)
    {
        idx += 1;
    } // while
    return idx < len && expression[idx] == delim ? idx : 0;
} // match_until


static size_t
match_variant(size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len, char const expression[static restrict len],
    GVA_Variant variant[static restrict 1])
{
    size_t idx = match_location(len, expression, &variant->start, &variant->end);
    if (idx == 0)
    {
        return 0;  // expected location
    } // if

    size_t tok = 0;
    if ((tok = match(len - idx, expression + idx, "dup")))
    {
        idx += tok;  // "dup"
        if (variant->end <= variant->start)
        {
            return 0;  // invalid range
        } // if

        tok = match_sequence(len - idx, expression + idx);
        if (tok)
        {
            if (tok != variant->end - variant->start)
            {
                return 0;  // inconsistent deleted length
            } // if
            if (memcmp(reference + variant->start, expression + idx, tok) != 0)
            {
                return 0;  // sequence not found in reference
            } // if
            idx += tok;
        } // if

        variant->sequence = (GVA_String) {variant->end - variant->start, reference + variant->start};
        return idx;
    } // if

    if ((tok = match(len - idx, expression + idx, "del")))
    {
        idx += tok;  // "del"
        if (variant->end <= variant->start)
        {
            return 0;  // invalid range
        } // if

        tok = match_sequence(len - idx, expression + idx);
        if (tok)
        {
            if (tok != variant->end - variant->start)
            {
                return 0;  // inconsistent deleted length
            } // if
            if (memcmp(reference + variant->start, expression + idx, tok) != 0)
            {
                return 0;  // sequence not found in reference
            } // if
            idx += tok;
        } // if

        if ((tok = match(len - idx, expression + idx, "ins")))
        {
            idx += tok;  // "ins"
            tok = match_sequence(len - idx, expression + idx);
            if (tok == 0)
            {
                return 0;  // missing insertion
            } // if
            variant->sequence = (GVA_String) {tok, expression + idx};
            idx += tok;
        } // if
        return idx;
    } // if

    if ((tok = match(len - idx, expression + idx, "ins")))
    {
        idx += tok;  // "ins"
        tok = match_sequence(len - idx, expression + idx);
        if (tok == 0)
        {
            return 0;  // missing insertion
        } // if
        variant->sequence = (GVA_String) {tok, expression + idx};
        idx += tok;
        return idx;
    } // if

    if ((tok = match_sequence(len - idx, expression + idx)))
    {
        if (tok != variant->end - variant->start)
        {
            return 0;  // inconsistent deleted length
        } // if
        if (memcmp(reference + variant->start, expression + idx, tok) != 0)
        {
            return 0;  // sequence not found in reference
        } // if
        idx += tok;
    } // if

    if (expression[idx] == '>')
    {
        idx += 1;  // '>'
        tok = match_sequence(len - idx, expression + idx);
        if (tok == 0)
        {
            return 0;  // missing insertion
        } // if
        variant->sequence = (GVA_String) {tok, expression + idx};
        idx += tok;
        return idx;
    } // if

    if (expression[idx] == '=')
    {
        idx += 1;  // '='
        return idx;
    } // if

    return 0;  // unexpected variant type
} // match_variant


GVA_Variant*
gva_parse_hgvs(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len, char const expression[static restrict len])
{
    // reference sequence identifier
    size_t idx = match_until(len, expression, ':');
    if (idx >= len)
    {
        return NULL;  // expected ':'
    } // if
    idx += 1;  // ':'

    // optional "g."
    idx += match(len - idx, expression + idx, "g.");
    if (idx >= len)
    {
        return NULL;  // expected variant
    } // if

    if (expression[idx] == '=')
    {
        return NULL;  // empty variant
    } // if

    if (expression[idx] == '[')
    {
        // TODO
        return NULL;
    } // if

    GVA_Variant variant = {0};
    size_t const tok = match_variant(len_ref, reference, len - idx, expression + idx, &variant);
    fprintf(stderr, "HIER %zu\n", tok);

    if (tok == 0 || idx + tok != len)
    {
        return NULL;  // expected variant
    } // if

    return NULL;
} // gva_parse_hgvs


size_t
gva_parse_spdi(size_t const len, char const expression[static restrict len],
    GVA_Variant variants[static restrict 1])
{
    variants[0] = (GVA_Variant) {0};
    // reference sequence identifier
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
