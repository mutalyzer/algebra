// NOT FREESTANDING
#ifndef GVA_UTILS_H
#define GVA_UTILS_H


#include <stddef.h>     // size_t
#include <stdio.h>      // FILE


#include "allocator.h"  // GVA_Allocator
#include "types.h"      // GVA_String


GVA_String
gva_fasta_sequence(GVA_Allocator const allocator, FILE* const restrict stream);


#endif // GVA_UTILS_H
