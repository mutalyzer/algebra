// NOT FREESTANDING
#include <stddef.h>     // NULL, size_t
#include <stdlib.h>     // rand
#include <stdio.h>      // FILE, fgets
#include <string.h>     // strcspn

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/string.h"      // GVA_String
#include "../include/utils.h"       // gva_fasta_sequence


GVA_String
gva_fasta_sequence(GVA_Allocator const allocator, FILE* const stream)
{
    static size_t const FASTA_LINE_SIZE = 65536;
    size_t capacity = FASTA_LINE_SIZE;
    GVA_String seq = {0, allocator.allocate(allocator.context, NULL, 0, capacity)};
    if (seq.str == NULL)
    {
        return seq;
    } // if

    while (fgets((char*) seq.str + seq.len, FASTA_LINE_SIZE, stream) != NULL)
    {
        if (seq.str[seq.len] == '>')
        {
            if (seq.len == 0)
            {
                continue;
            } // if
            break;
        } // if
        seq.len += strcspn(seq.str + seq.len, "\n");

        if (capacity - FASTA_LINE_SIZE < seq.len)
        {
            seq.str = allocator.allocate(allocator.context, (char*) seq.str, capacity, capacity + FASTA_LINE_SIZE);
            if (seq.str == NULL)
            {
                seq.len = 0;
                return seq;
            } // if
            capacity += FASTA_LINE_SIZE;
        } // if
    } // while
    seq.str = allocator.allocate(allocator.context, (char*) seq.str, capacity, seq.len);
    return seq;
} // gva_fasta_sequence


static inline void
random_sequence(size_t const len, char sequence[static len])
{
    for (size_t i = 0; i < len; ++i)
    {
        sequence[i] = "AC"[rand() % 2];
    } // for
} // random_sequence


GVA_String
gva_random_sequence(GVA_Allocator const allocator, size_t const min_length, size_t const max_length)
{
    if (max_length < min_length)
    {
        return (GVA_String) {0, NULL};
    } // if

    size_t const len = rand() % (max_length - min_length + 1) + min_length;
    GVA_String sequence = {len, allocator.allocate(allocator.context, NULL, 0, len)};
    if (sequence.str == NULL)
    {
        return (GVA_String) {0, NULL};
    } // if
    random_sequence(len, (char*) sequence.str);
    return sequence;
} // gva_random_sequence
