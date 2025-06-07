// NOT FREESTANDING
#ifndef GVA_UTILS_H
#define GVA_UTILS_H


#include <stddef.h>     // size_t
#include <stdio.h>      // FILE


#include "allocator.h"  // GVA_Allocator


#define GVA_VARIANT_PRINT(variant, observed) variant.start, variant.end, (int) variant.obs_end - variant.obs_start, observed + variant.obs_start
#define GVA_VARIANT_FMT "%u:%u/%.*s"


typedef struct
{
    size_t len;
    char* str;
} GVA_String;


GVA_String
gva_fasta_sequence(GVA_Allocator const allocator, FILE* const restrict stream);


#endif // GVA_UTILS_H