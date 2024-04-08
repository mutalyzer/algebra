#include <stdbool.h>    // bool
#include <stddef.h>     // size_t
#include <stdint.h>     // uint32_t

#include "../include/parser.h"      // va_parse_*
#include "../include/variant.h"     // VA_Variant


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
match_number(size_t const len, char const expression[static len],
             uint32_t* const num)
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
match_location(size_t const len, char const expression[static len],
               uint32_t* const start, uint32_t* const end)
{
    size_t idx = match_number(len, expression, start);
    if (idx == 0)
    {
        return 0;  // expected number
    } // if
    if (idx < len && expression[idx] == '_')
    {
        idx += 1;
        size_t const tok = match_number(len - idx, expression + idx, end);
        if (tok == 0)
        {
            return 0;  // expected number
        } // if
        if (*start >= *end)
        {
            return 0;  // invalid range
        } // if
        return idx + tok;
    } // if
    *end = *start;
    return idx;
} // match_location


static size_t
match(size_t const len, char const expression[static len], char const word[static 1])
{
    size_t idx = 0;
    while (idx < len && word[idx] != '\0' && expression[idx] == word[idx])
    {
        idx += 1;
    } // while
    return word[idx] == '\0' ? idx : 0;
} // match


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


static size_t
match_variant(size_t const len, char const expression[static len],
              VA_Variant* const variant)
{
    size_t idx = match_location(len, expression, &variant->start, &variant->end);
    if (idx == 0)
    {
        return 0;  // expected location
    } // if

    size_t tok = match(len - idx, expression + idx, "del");
    if (tok != 0)
    {
        if (variant->start == 0)
        {
            return 0;  // invalid location
        } // if
        idx += tok;
        variant->start -= 1;

        // optional sequence
        idx += match_sequence(len - idx, expression + idx);

        tok = match(len - idx, expression + idx, "ins");
        if (tok != 0)
        {
            idx += tok;

            variant->obs_start = idx;  // OVERFLOW
            tok = match_sequence(len - idx, expression + idx);
            if (tok == 0)
            {
                return 0;  // expected sequence;
            } // if
            idx += tok;

            variant->obs_end = idx;  // OVERFLOW
        } // if

        return idx;
    } // if

    tok = match(len - idx, expression + idx, "ins");
    if (tok != 0)
    {
        if (variant->end - variant->start != 1)
        {
            return 0;  // expected 1-interval
        } // if
        idx += tok;
        variant->end = variant->start;

        variant->obs_start = idx;  // OVERFLOW
        tok = match_sequence(len - idx, expression + idx);
        if (tok == 0)
        {
            return 0;  // expected sequence
        } // if
        idx += tok;
        variant->obs_end = idx;  // OVERFLOW

        return idx;
    } // if

    if (variant->start == 0)
    {
        return 0;  // invalid location
    } // if
    variant->start -= 1;

    // optional sequence
    idx += match_sequence(len - idx, expression + idx);

    if (idx < len && expression[idx] == '>')
    {
        idx += 1;  // '>'

        variant->obs_start = idx;  // OVERFLOW
        tok = match_sequence(len - idx, expression + idx);
        if (tok == 0)
        {
            return 0;  // expected sequence
        } // if
        idx += tok;
        variant->obs_end = idx;  // OVERFLOW

        return idx;
    } // if

    return 0;
} // match_variant


size_t
va_parse_hgvs(size_t const len, char const expression[static len],
              size_t const len_var, VA_Variant variants[static len_var])
{
    size_t idx = match_until(len, expression, ':');
    if (idx >= len)
    {
        return 0;  // expected ':'
    } // if
    if (idx > 0)
    {
        idx += 1;  // ':'
    } // if
    idx += match(len - idx, expression + idx, "g.");

    if (idx < len && expression[idx] == '[')
    {
        idx += 1;

        size_t tok = match_variant(len - idx, expression + idx, &variants[0]);
        if (tok == 0)
        {
            return 0;  // expected variant
        } // if
        idx += tok;

        size_t count = 1;
        while (idx < len && expression[idx] == ';')
        {
            idx += 1;  // ';'

            if (count >= len_var)
            {
                return 0;  // too many variants
            } // if

            tok = match_variant(len - idx, expression + idx, &variants[count]);
            if (tok == 0)
            {
                return 0;  // expected variant
            } // if
            idx += tok;
            count += 1;
        } // while

        if (idx >= len || expression[idx] != ']')
        {
            return 0;  // expected ']'
        } // if
        idx += 1;  // ']'

        return idx == len ? count : 0;
    } // if

    if (idx < len && expression[idx] == '=')
    {
        idx += 1;
        return idx == len;
    } // if

    size_t tok = match_variant(len - idx, expression + idx, &variants[0]);
    if (tok == 0)
    {
        return 0;  // expected variant
    } // if
    idx += tok;

    return idx == len;
} // va_parse_hgvs


size_t
va_parse_spdi(size_t const len, char const expression[static len],
              VA_Variant variants[static 1])
{
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
        if (tok == 0)
        {
            return 0;  // expected number (length) or sequence
        } // if
        variants[0].end = tok;  // OVERFLOW
    } // if
    idx += tok;
    variants[0].end += variants[0].start;

    if (idx >= len || expression[idx] != ':')
    {
        return 0;  // expected ':'
    } // if
    idx += 1;

    variants[0].obs_start = idx;  // OVERFLOW
    tok = match_sequence(len - idx, expression + idx);
    if (tok == 0)
    {
        return 0;  // expected sequence
    } // if
    idx += tok;
    variants[0].obs_end = idx;

    return idx == len;
} // va_parse_spdi
