#ifndef GVA_VARIANT_H
#define GVA_VARIANT_H


#include <stddef.h>     // size_t


#include "allocator.h"  // GVA_Allocator
#include "types.h"      // gva_uint, GVA_String


#define GVA_VARIANT_PRINT(variant) variant.start, variant.end, (int) variant.sequence.len, variant.sequence.str
#define GVA_VARIANT_FMT "%u:%u/%.*s"


typedef struct
{
    gva_uint   start;
    gva_uint   end;
    GVA_String sequence;
} GVA_Variant;


size_t
gva_parse_spdi(size_t const len, char const expression[static restrict len],
    GVA_Variant variants[static restrict 1]);


size_t
gva_variant_length(GVA_Variant const variant);


GVA_String
gva_patch(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const n, GVA_Variant const variants[static restrict n]);


#endif // GVA_VARIANT_H
