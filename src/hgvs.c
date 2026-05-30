// NOT FREESTANDING
#include <errno.h>      // errno
#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // ...
#include <stdlib.h>     // ...
#include <string.h>     // ...

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/parser.h"      // gva_parse_*
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/string.h"      // GVA_STRING_*
#include "../include/utils.h"       // gva_fasta_sequence_blob
#include "../include/variant.h"     // GVA_Variant

#include "array.h"  // ARRAY_DESTROY


#define LINE_SIZE 4096


static inline size_t
repeats(GVA_Allocator const allocator,
    size_t const len, char const word[static len])
{
    size_t* lps = allocator.allocate(allocator.context, NULL, 0, len * sizeof(*lps));
    if (lps == NULL)
    {
        return 0;  // OOM
    } // if

    lps[0] = 0;
    size_t length = 0;
    size_t idx = 1;
    while (idx < len)
    {
        if (word[idx] == word[length])
        {
            length += 1;
            lps[idx] = length;
            idx += 1;
        } // if
        else if (length > 0)
        {
            length = lps[length - 1];
        } // if
        else
        {
            lps[idx] = 0;
            idx += 1;
        } // else
    } // while

    lps = allocator.allocate(allocator.context, lps, len * sizeof(*lps), 0);
    return len - length;
} // repeats


static inline size_t
hgvs_location(size_t const len, char buffer[static len],
    size_t const start, size_t const end)
{
    if (end - start == 0)
    {
        return snprintf(buffer, len, "%zu_%zu", start, start + 1);
    } // if
    if (end - start == 1)
    {
        return snprintf(buffer, len, "%zu", start + 1);
    } // if
    return snprintf(buffer, len, "%zu_%zu", start + 1, end);
} // hgvs_location


static size_t
hgvs_variant(size_t const len, char buffer[static len],
    size_t const len_ref, char const reference[static restrict len_ref],
    GVA_Variant const variant)
{
    size_t const inserted = repeats(gva_std_allocator, variant.sequence.len, variant.sequence.str);
    size_t const deleted = repeats(gva_std_allocator, variant.end - variant.start, reference + variant.start);

    fprintf(stderr, "%zu vs %zu\n", inserted, deleted);


    size_t idx = hgvs_location(len, buffer, variant.start, variant.end);

    if (variant.end - variant.start == 0)
    {
        if (variant.sequence.len == 0)
        {
            return idx + snprintf(buffer + idx, len - idx, "=");
        } // if
        return idx + snprintf(buffer + idx, len - idx, "ins" GVA_STRING_FMT, GVA_STRING_PRINT(variant.sequence));
    } // if

    if (variant.end - variant.start == 1)
    {
        if (variant.sequence.len == 0)
        {
            return idx + snprintf(buffer + idx, len - idx, "del");
        } // if
        if (variant.sequence.len == 1)
        {
            return idx + snprintf(buffer + idx, len - idx, "%c>" GVA_STRING_FMT, reference[variant.start], GVA_STRING_PRINT(variant.sequence));
        } // if
        return idx + snprintf(buffer + idx, len - idx, "delins" GVA_STRING_FMT, GVA_STRING_PRINT(variant.sequence));
    } // if

    if (variant.sequence.len == 0)
    {
        return idx + snprintf(buffer + idx, len - idx, "del");
    } // if

    return idx + snprintf(buffer + idx, len - idx, "delins" GVA_STRING_FMT, GVA_STRING_PRINT(variant.sequence));
} // hgvs_variant


static size_t
to_hgvs(size_t const len, char buffer[static len],
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const n, GVA_Variant const variants[static restrict n])
{
    if (n == 0)
    {
        return snprintf(buffer, len, "=");
    } // if

    if (n == 1)
    {
        return hgvs_variant(len, buffer, len_ref, reference, variants[0]);
    } // if

    size_t idx = snprintf(buffer, len, "[");
    for (size_t i = 0; i < n; ++i)
    {
        idx += hgvs_variant(len - idx, buffer + idx, len_ref, reference, variants[i]);
        if (i < n - 1)
        {
            idx += snprintf(buffer + idx, len - idx, ";");
        } // if
    } // for
    return idx + snprintf(buffer + idx, len - idx, "]");
} // to_hgvs


int
hgvs_main(int argc, char* argv[static argc])
{
    static char buffer[1024] = {'\0'};
    GVA_String const ref = {10, "AAATAATATAATAATTTAT"};
    GVA_Variant const variant = {2, 13, {6, "ATAATA"}};
    fprintf(stderr, "%zu\n", to_hgvs(1024, buffer, ref.len, ref.str, 1, &variant));

    fprintf(stderr, "%s\n", buffer);

    return EXIT_SUCCESS;


    if (argc < 2)
    {
        fprintf(stderr, "usage %s reference.blob\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    errno = 0;
    FILE* stream = fopen(argv[1], "r");
    if (stream == NULL)
    {
        perror("fopen()");
        return EXIT_FAILURE;
    } // if

    GVA_String reference = {0};
    reference = gva_fasta_sequence_blob(gva_std_allocator, stream);
    fclose(stream);

    fprintf(stderr, "reference length: %zu\n", reference.len);

    size_t line_count = 0;
    static char line[LINE_SIZE] = {0};
    while (fgets(line, sizeof(line), stdin) != NULL)
    {
        line_count += 1;

        size_t idx = 0;
        size_t tok = 0;

        tok = strcspn(line + idx, "\t ");
        idx += tok + 1;  // ignore identifier

        tok = strcspn(line + idx, "\t ");
        GVA_Variant variant = {0};
        if (gva_parse_spdi(tok, line + idx, &variant) == 0)
        {
            fprintf(stderr, "%zu: ERROR SPDI: %s", line_count, line);
            continue;
        } // if
        idx += tok + 1;

        tok = strcspn(line + idx, "\n\t ");
        GVA_HGVS_Allele allele = gva_parse_hgvs(gva_std_allocator, reference.len, reference.str, tok, line + idx );
        if (!allele.interpretable)
        {
            fprintf(stderr, "%zu: ERROR HGVS: %s", line_count, line);
            continue;
        } // if
        idx += tok + 1;

        ARRAY_DESTROY(gva_std_allocator, allele.inserted);
        ARRAY_DESTROY(gva_std_allocator, allele.variants);
    } //while

    gva_string_destroy(gva_std_allocator, reference);

    return EXIT_SUCCESS;
} // hgvs_main
