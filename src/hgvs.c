// NOT FREESTANDING
#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // snprintf
#include <string.h>     // memcmp

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/string.h"      // GVA_STRING_*
#include "../include/types.h"       // GVA_UINT_FMT
#include "../include/variant.h"     // GVA_Variant


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


static size_t
hgvs_variant(size_t const len, char buffer[static len],
    size_t const len_ref, char const reference[static restrict len_ref],
    GVA_Variant const variant)
{
    fprintf(stderr, GVA_STRING_FMT "\n", GVA_STRING_PRINT(((GVA_String) {variant.end - variant.start, reference + variant.start})));
    fprintf(stderr, GVA_VARIANT_FMT "\n", GVA_VARIANT_PRINT(variant));

    size_t len_inserted = repeats(gva_std_allocator, variant.sequence.len, variant.sequence.str);
    if (len_inserted > 0)
    {
        size_t cnt_inserted = variant.sequence.len / len_inserted;
        size_t rem_inserted = variant.sequence.len % len_inserted;

        size_t len_deleted = repeats(gva_std_allocator, variant.end - variant.start, reference + variant.start);
        if (len_deleted > 0)
        {
            size_t cnt_deleted = (variant.end - variant.start) / len_deleted;
            size_t rem_deleted = (variant.end - variant.start) % len_deleted;

            char const* inserted = variant.sequence.str;
            char const* deleted = reference + variant.start;

            fprintf(stderr, "%zu: " GVA_STRING_FMT " %zu %zu\n", len_inserted,
                GVA_STRING_PRINT(((GVA_String) {len_inserted, inserted})),
                cnt_inserted, rem_inserted);
            fprintf(stderr, "%zu: " GVA_STRING_FMT " %zu %zu\n", len_deleted,
                GVA_STRING_PRINT(((GVA_String) {len_deleted, deleted})),
                cnt_deleted, rem_deleted);

            if (len_inserted < len_deleted &&
                memcmp(deleted, inserted, len_inserted) == 0)
            {
                fprintf(stderr, "if\n");
                len_inserted = len_deleted;
                cnt_inserted = 1;
                rem_inserted = rem_deleted;
            } // if
            else if (len_deleted < len_inserted &&
                memcmp(deleted, inserted, len_deleted) == 0)
            {
                fprintf(stderr, "else if\n");
                len_deleted = len_inserted;
                cnt_deleted = 1;
                rem_deleted = rem_inserted;
            } // if

            if (memcmp(deleted, inserted, len_deleted) == 0)
            {
                // duplictation
                if (cnt_deleted == 1 && cnt_inserted == 2)
                {
                    if (variant.end - variant.start == 0)
                    {
                        return snprintf(buffer, len, "%zudup",
                            variant.start + rem_deleted + 1);
                    } // if
                    return snprintf(buffer, len, "%zu_" GVA_UINT_FMT "dup",
                        variant.start + rem_deleted + 1, variant.end);
                } // if

                // repeat
                return snprintf(buffer, len, "%zu_" GVA_UINT_FMT GVA_STRING_FMT "[%zu]",
                    variant.start + rem_deleted + 1, variant.end,
                    GVA_STRING_PRINT(((GVA_String) {len_inserted, inserted + rem_inserted})),
                    cnt_inserted);
            } // if
        } // if
    } // if








    if (variant.end - variant.start == 0)
    {
        size_t const idx = snprintf(buffer, len, GVA_UINT_FMT "_" GVA_UINT_FMT, variant.start, variant.start + 1);
        if (variant.sequence.len == 0)
        {
            return idx + snprintf(buffer + idx, len - idx, "=");
        } // if
        return idx + snprintf(buffer + idx, len - idx, "ins" GVA_STRING_FMT, GVA_STRING_PRINT(variant.sequence));
    } // if

    if (variant.end - variant.start == 1)
    {
        size_t const idx = snprintf(buffer, len, GVA_UINT_FMT, variant.start + 1);
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

    size_t const idx = snprintf(buffer, len, GVA_UINT_FMT "_" GVA_UINT_FMT, variant.start + 1, variant.end);
    if (variant.sequence.len == 0)
    {
        return idx + snprintf(buffer + idx, len - idx, "del");
    } // if

    return idx + snprintf(buffer + idx, len - idx, "delins" GVA_STRING_FMT, GVA_STRING_PRINT(variant.sequence));
} // hgvs_variant


size_t
gva_to_hgvs(size_t const len, char buffer[static restrict len],
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
} // gva_to_hgvs
