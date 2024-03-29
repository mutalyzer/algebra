#include <stdbool.h>    // bool
#include <stddef.h>     // size_t
#include <stdint.h>     // uint32_t

#include "../include/parser.h"  // va_parse_*


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
    if (expression[idx] == '_')
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


size_t
va_parse_hgvs(size_t const len, char const expression[static len])
{
    return match(len, expression, "hgvs");
} // va_parse_hgvs


size_t
va_parse_spdi(size_t const len, char const expression[static len])
{
    return match(len, expression, "spdi");
} // va_parse_spdi
