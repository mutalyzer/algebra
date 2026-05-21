#include <stdbool.h>    // bool, true
#include <stddef.h>     // NULL, size_t
#include <string.h>     // memcmp, mempcy

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/parser.h"      // GVA_HGVS_Allele, gva_parse_*
#include "../include/string.h"      // GVA_String
#include "../include/types.h"       // gva_uint
#include "../include/variant.h"     // GVA_Variant
#include "array.h"      // ARRAY_*, array_*


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


inline char
gva_complement(char const ch)
{
    if (ch == 'A')
    {
        return 'T';
    } // if
    if (ch == 'C')
    {
        return 'G';
    } // if
    if (ch == 'G')
    {
        return 'C';
    } // if
    if (ch == 'T')
    {
        return 'A';
    } // if
    return ch;
} // gva_complement


static inline size_t
match(size_t const len, char const expression[static restrict len],
    size_t const len_word, char const word[static restrict len_word])
{
    size_t idx = 0;
    while (idx < len && idx < len_word && expression[idx] == word[idx])
    {
        idx += 1;
    } // while
    return idx == len_word ? idx : 0;
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


static inline size_t
match_location(size_t const len, char const expression[static restrict len],
    gva_uint start[static restrict 1], gva_uint end[static restrict 1])
{
    size_t idx = 0;
    size_t tok = 0;

    if (!(tok = match_number(len - idx, expression + idx, start)))
    {
        //fprintf(stderr, "expected number\n");
        return 0;
    } // if
    idx += tok;

    *end = *start;
    if (*start > 0)
    {
        *start -= 1;
    } // if
    if (idx < len && expression[idx] == '_')
    {
        idx += 1;  // '_'
        *end = 0;
        if (!(tok = match_number(len - idx, expression + idx, end)))
        {
            //fprintf(stderr, "expected number\n");
            return 0;
        } // if
        if (*end < *start)
        {
            //fprintf(stderr, "invalid range\n");
            return 0;
        } // if
        idx += tok;
    } // if
    return idx;
} // match_location


static size_t
match_inserted_part(GVA_Allocator const allocator,
    GVA_HGVS_Allele allele[static restrict 1],
    size_t const len, char const expression[static restrict len])
{
    size_t idx = 0;
    size_t tok = 0;

    if (!(tok = match_sequence(len - idx, expression + idx)))
    {
        //fprintf(stderr, "expected sequence\n");
        return 0;  // expected sequence
    } // if
    GVA_String const sequence = {tok, expression + idx};
    idx += tok;

    if (idx < len && expression[idx] == '[')
    {
        idx += 1;  // '['
        gva_uint count = 0;

        if (!(tok = match_number(len - idx, expression + idx, &count)))
        {
            //fprintf(stderr, "expected number\n");
            return 0;  // expected number
        } // if
        idx += tok;

        if (idx < len && expression[idx] == ']')
        {
            idx += 1;  // ']'

            if (count > 0)
            {
                allele->inserted = array_ensure(allocator, allele->inserted, 1, sequence.len * count);  // OVERFLOW
                for (gva_uint i = 0; i < count; ++i)
                {
                    memcpy(allele->inserted + array_length(allele->inserted) + i * sequence.len, sequence.str, sequence.len);
                } // for
                array_header(allele->inserted)->length += sequence.len * count;
            } // if
            return idx;
        } // if

        //fprintf(stderr, "expected ']'\n");
        return 0;  // expected ']'
    } // if

    allele->inserted = array_ensure(allocator, allele->inserted, 1, sequence.len);
    memcpy(allele->inserted + array_length(allele->inserted), sequence.str, sequence.len);
    array_header(allele->inserted)->length += sequence.len;

    return idx;
} // match_inserted_part


static size_t
match_insertion(GVA_Allocator const allocator,
    GVA_HGVS_Allele allele[static restrict 1],
    size_t const len, char const expression[static restrict len])
{
    size_t idx = 0;
    size_t tok = 0;

    if (idx < len && expression[idx] == '[')
    {
        idx += 1;  // '['
        while (idx < len)
        {
            if (!(tok = match_inserted_part(allocator, allele, len - idx, expression + idx)))
            {
                //fprintf(stderr, "expected inserted part\n");
                return 0;  // expected inserted part
            } // if
            idx += tok;

            if (idx < len && expression[idx] != ';')
            {
                break;
            } // if
            idx += 1;  // ';'
        } // while

        if (idx < len && expression[idx] == ']')
        {
            idx += 1;  // ']'
            return idx;
        } // if

        //fprintf(stderr, "expected ']'\n");
        return 0;  // expected ']'
    } // if

    if (!(tok = match_inserted_part(allocator, allele, len - idx, expression + idx)))
    {
        //fprintf(stderr, "expected inserted part\n");
        return 0;  // expected inserted part
    } // if
    idx += tok;
    return idx;
} // match_insertion


static size_t
match_variant(GVA_Allocator const allocator,
    GVA_HGVS_Allele allele[static restrict 1],
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len, char const expression[static restrict len],
    GVA_Variant variant[static restrict 1])
{
    size_t idx = 0;
    size_t tok = 0;

    if (!(tok = match_location(len - idx, expression + idx, &variant->start, &variant->end)))
    {
        //fprintf(stderr, "expected location\n");
        return 0;  // expected location
    } // if

    if (variant->end >= len_ref)
    {
        //fprintf(stderr, "invalid location in reference\n");
        return 0;  // invalid location in reference
    } // if
    idx += tok;

    // deletion and optionally insertion
    if ((tok = match(len - idx, expression + idx, 3, "del")))
    {
        if (variant->start == variant->end)
        {
            //fprintf(stderr, "invalid deleted range\n");
            return 0;  // invalid deleted range
        } // if
        idx += tok;

        // optional number
        gva_uint count = 0;
        if ((tok = match_number(len - idx, expression + idx, &count)))
        {
            idx += tok;
            if (count != variant->end - variant->start)
            {
                //fprintf(stderr, "inconsistent deleted length\n");
                return 0;  // inconsistent deleted length
            } // if
        } // if
        // optional sequence
        else if ((tok = match_sequence(len - idx, expression + idx)))
        {
            if (variant->end - variant->start != tok)
            {
                //fprintf(stderr, "inconsistent deleted length\n");
                return 0;  // inconsistent deleted length
            } // if
            if (memcmp(reference + variant->start, expression + idx, tok))
            {
                //fprintf(stderr, "sequence not found in reference\n");
                return 0;  // sequence not found in reference
            } // if
            idx += tok;
        } // if

        // optional insertion
        if ((tok = match(len - idx, expression + idx, 3, "ins")))
        {
            idx += tok;
            size_t const start = array_length(allele->inserted);
            if ((tok = match_insertion(allocator, allele, len - idx, expression + idx)))
            {
                idx += tok;
                variant->sequence = (GVA_String) {array_length(allele->inserted) - start, (void*) start};
            } // if
        } // if

        return idx;
    } // if

    // duplication
    if ((tok = match(len - idx, expression + idx, 3, "dup")))
    {
        if (variant->start == variant->end)
        {
            //fprintf(stderr, "invalid duplicated range\n");
            return 0;  // invalid duplicated range
        } // if
        idx += tok;

        // optional sequence
        if ((tok = match_sequence(len - idx, expression + idx)))
        {
            if (variant->end - variant->start != tok)
            {
                //fprintf(stderr, "inconsistent duplicated length\n");
                return 0;  // inconsistent duplicated length
            } // if
            if (memcmp(reference + variant->start, expression + idx, tok))
            {
                //fprintf(stderr, "sequence not found in reference\n");
                return 0;  // sequence not found in reference
            } // if
            idx += tok;
        } // if

        size_t const length = variant->end - variant->start;
        size_t const start = array_length(allele->inserted);
        allele->inserted = array_ensure(allocator, allele->inserted, 1, length);
        memcpy(allele->inserted + start, reference + variant->start, length);
        variant->sequence = (GVA_String) {length, (void*) start};
        return idx;
    } // if

    // inversion
    if ((tok = match(len - idx, expression + idx, 3, "inv")))
    {
        if (variant->start == variant->end)
        {
            //fprintf(stderr, "invalid inverted range\n");
            return 0;  // invalid inverted range
        } // if
        idx += tok;

        // optional sequence
        if ((tok = match_sequence(len - idx, expression + idx)))
        {
            if (variant->end - variant->start != tok)
            {
                //fprintf(stderr, "inconsistent inverted length\n");
                return 0;  // inconsistent inverted length
            } // if

            for (size_t i = 0; i < tok; ++i)
            {
                if (reference[variant->end - i - 1] != gva_complement(expression[idx + i]))
                {
                    //fprintf(stderr, "complement not found in reference\n");
                    return 0;  // complement not found in reference
                } // if
            } // for

            idx += tok;
        } // if

        size_t const length = variant->end - variant->start;
        size_t const start = array_length(allele->inserted);
        allele->inserted = array_ensure(allocator, allele->inserted, 1, length);
        for (size_t i = 0; i < length; ++i)
        {
            allele->inserted[start + i] = gva_complement(reference[variant->end - i - 1]);
        } // for
        variant->sequence = (GVA_String) {length, (void*) start};
        return idx;
    } // if

    // insertion
    if ((tok = match(len - idx, expression + idx, 3, "ins")))
    {
        if (variant->end - variant->start > 2)
        {
            //fprintf(stderr, "invalid range for insertion\n");
            return 0;  // invalid range for insertion
        } // if
        idx += tok;
        variant->end = variant->start;

        size_t const start = array_length(allele->inserted);
        if ((tok = match_insertion(allocator, allele, len - idx, expression + idx)))
        {
            idx += tok;
            variant->sequence = (GVA_String) {array_length(allele->inserted) - start, (void*) start};
        } // if
        return idx;
    } // if

    GVA_String sequence = {0};
    // optional sequence
    if ((tok = match_sequence(len - idx, expression + idx)))
    {
        sequence = (GVA_String) {tok, expression + idx};
        idx += tok;
    } // if

    // substitution
    if (idx < len && expression[idx] == '>')
    {
        idx += 1;  // '>'

        if (sequence.len > 0)
        {
            if (sequence.len != variant->end - variant->start)
            {
                //fprintf(stderr, "inconsistent sequence length\n");
                return 0;  // inconsistent sequence length
            } // if
            if (memcmp(reference + variant->start, sequence.str, sequence.len))
            {
                //fprintf(stderr, "sequence not found in reference\n");
                return 0;  // sequence not found in reference
            } // if
        } // if

        size_t const start = array_length(allele->inserted);
        if (!(tok = match_inserted_part(allocator, allele, len - idx, expression + idx)))
        {
            //fprintf(stderr, "expected inserted part\n");
            return 0;  // expected inserted part
        } // if
        idx += tok;
        variant->sequence = (GVA_String) {array_length(allele->inserted) - start, (void*) start};
        return idx;
    } // if

    // empty variant
    if (idx < len && expression[idx] == '=')
    {
        idx += 1;  // '='
        size_t const length = variant->end - variant->start;
        size_t const start = array_length(allele->inserted);
        allele->inserted = array_ensure(allocator, allele->inserted, 1, length);
        memcpy(allele->inserted + start, reference + variant->start, length);
        variant->sequence = (GVA_String) {length, (void*) start};
        return idx;
    } // if

    // repeats
    if (idx < len && expression[idx] == '[')
    {
        idx += 1;  // '['

        variant->end -= 1;
        if (sequence.len == 0)
        {
            //fprintf(stderr, "expected sequence\n");
            return 0;  // expected sequence
        } // if

        gva_uint count = 0;
        if (!(tok = match_number(len - idx, expression + idx, &count)))
        {
            //fprintf(stderr, "expected number\n");
            return 0;  // expected number
        } // if
        idx += tok;

        if (idx < len && expression[idx] == ']')
        {
            idx += 1;  // ']'
            // NCBI style repeat
            if (variant->end - variant->start == 0)
            {
                size_t found = 0;
                while (match(len_ref - variant->start + found * sequence.len, reference + variant->start + found * sequence.len, sequence.len, sequence.str))
                {
                    found += 1;
                } // while
                if (found == 0)
                {
                    //fprintf(stderr, "sequence not found in reference\n");
                    return 0;  // sequence not found in reference
                } // if
                variant->end = variant->start + found * sequence.len;

                size_t const start = array_length(allele->inserted);
                if (count > 0)
                {
                    allele->inserted = array_ensure(allocator, allele->inserted, 1, sequence.len * count);  // OVERFLOW
                    for (gva_uint i = 0; i < count; ++i)
                    {
                        memcpy(allele->inserted + array_length(allele->inserted) + i * sequence.len, sequence.str, sequence.len);
                    } // for
                    array_header(allele->inserted)->length += sequence.len * count;
                } // if
                variant->sequence = (GVA_String) {sequence.len * count, (void*) start};
                return idx;
            } // if

            size_t const start = array_length(allele->inserted);
            if (count > 0)
            {
                allele->inserted = array_ensure(allocator, allele->inserted, 1, sequence.len * count);  // OVERFLOW
                for (gva_uint i = 0; i < count; ++i)
                {
                    memcpy(allele->inserted + array_length(allele->inserted) + i * sequence.len, sequence.str, sequence.len);
                } // for
                array_header(allele->inserted)->length += sequence.len * count;
            } // if
            variant->sequence = (GVA_String) {sequence.len * count, (void*) start};
            return idx;
        } // if

        //fprintf(stderr, "expected ']'\n");
        return 0;  // expected ']'
    } // if

    //fprintf(stderr, "unsupported variant\n");
    return 0;  // unsupported variant
} // match_variant


GVA_HGVS_Allele
gva_parse_hgvs(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len, char const expression[static restrict len])
{
    size_t idx = 0;
    size_t tok = 0;

    // optional reference sequence identifier
    if ((tok = match_until(len - idx, expression + idx, ':')))
    {
        idx += tok + 1;  // ':';
    } // if

    // optional "g."
    idx += match(len - idx, expression + idx, 2, "g.");

    // empty variant
    if (idx < len && expression[idx] == '=')
    {
        idx += 1;  // '='
        if (idx == len)
        {
            return (GVA_HGVS_Allele) {.interpretable = true};
        } // if

        //fprintf(stderr, "expected end of expression\n");
        return (GVA_HGVS_Allele) {NULL};
    } // if

    GVA_HGVS_Allele allele = {NULL};

    // allele
    if (idx < len && expression[idx] == '[')
    {
        idx += 1;  // '['

        while (idx < len)
        {
            GVA_Variant variant = {0};
            if ((tok = match_variant(allocator, &allele, len_ref, reference, len - idx, expression + idx, &variant)))
            {
                idx += tok;
                ARRAY_APPEND(allocator, allele.variants, variant);
            } // if

            if (idx < len && expression[idx] != ';')
            {
                break;
            } // if
            idx += 1;  // ';'
        } // while

        if (idx < len && expression[idx] == ']')
        {
            idx += 1;  // ']'
            if (idx == len)
            {
                for (size_t i = 0; i < array_length(allele.variants); ++i)
                {
                    allele.variants[i].sequence.str = allele.inserted + (size_t) allele.variants[i].sequence.str;
                } // for
                allele.interpretable = true;
                return allele;
            } // if

            allele.variants = ARRAY_DESTROY(allocator, allele.variants);
            allele.inserted = ARRAY_DESTROY(allocator, allele.inserted);
            //fprintf(stderr, "expected end of expression\n");
            return (GVA_HGVS_Allele) {NULL};  // expected end of expression
        } // if

        allele.variants = ARRAY_DESTROY(allocator, allele.variants);
        allele.inserted = ARRAY_DESTROY(allocator, allele.inserted);
        //fprintf(stderr, "expected ']'\n");
        return (GVA_HGVS_Allele) {NULL};  // expected ']'
    } // if

    // single variant
    GVA_Variant variant = {0};
    if ((tok = match_variant(allocator, &allele, len_ref, reference, len - idx, expression + idx, &variant)))
    {
        idx += tok;
        if (idx == len)
        {
            variant.sequence.str = allele.inserted + (size_t) variant.sequence.str;
            ARRAY_APPEND(allocator, allele.variants, variant);
            allele.interpretable = true;
            return allele;
        } // if

        allele.inserted = ARRAY_DESTROY(allocator, allele.inserted);
        //fprintf(stderr, "expected end of expression\n");
        return (GVA_HGVS_Allele) {NULL};
    } // if

    //fprintf(stderr, "expected variant\n");
    return (GVA_HGVS_Allele) {NULL};  // expected variant
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
