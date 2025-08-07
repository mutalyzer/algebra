// NOT FREESTANDING
#ifndef GVA_UTILS_H
#define GVA_UTILS_H


#include <stddef.h>     // size_t
#include <stdio.h>      // FILE

#include "allocator.h"  // GVA_Allocator
#include "string.h"     // GVA_String


// Reads the first sequence from a stream in FASTA format. A single header
// line (starting with '>') is expected after which the sequence is taken
// as is only omitting '\n' characters until EOF or the next '>' at the
// beginning of a line.
// The result must be deallocated by the caller (e.g. `gva_string_destroy`).
GVA_String
gva_fasta_sequence(GVA_Allocator const allocator, FILE* const stream);


GVA_String
gva_random_sequence(GVA_Allocator const allocator,
    size_t const min_length, size_t const max_length);


#endif // GVA_UTILS_H
