// NOT FREESTANDING
#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // FILE, fgets
#include <string.h>     // strcspn


#include "../include/allocator.h"   // GVA_Allocator
#include "../include/utils.h"       // GVA_String, gva_fasta_sequence


GVA_String
gva_fasta_sequence(GVA_Allocator const allocator, FILE* const restrict stream)
{
    static size_t const FASTA_LINE_SIZE = 65536;
    size_t capacity = FASTA_LINE_SIZE;
    GVA_String seq = {allocator.allocate(allocator.context, NULL, 0, capacity), 0};
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
