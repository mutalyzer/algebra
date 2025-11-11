#ifndef GVA_VARIANT_H
#define GVA_VARIANT_H


#include <stdbool.h>    // bool
#include <stddef.h>     // size_t

#include "allocator.h"  // GVA_Allocator
#include "string.h"     // GVA_STRING_*, GVA_String
#include "types.h"      // GVA_UINT_FMT, gva_uint


#define GVA_VARIANT_PRINT(variant) variant.start, variant.end, GVA_STRING_PRINT(variant.sequence)
#define GVA_VARIANT_FMT GVA_UINT_FMT ":" GVA_UINT_FMT "/" GVA_STRING_FMT

#define GVA_VARIANT_PRINT_SPDI(ref_id, variant) ref_id, variant.start, (variant.end - variant.start), GVA_STRING_PRINT(variant.sequence)
#define GVA_VARIANT_FMT_SPDI "%s:" GVA_UINT_FMT ":" GVA_UINT_FMT ":" GVA_STRING_FMT


typedef struct
{
    gva_uint   start;
    gva_uint   end;
    GVA_String sequence;
} GVA_Variant;


size_t
gva_parse_spdi(size_t const len, char const expression[static restrict len],
    GVA_Variant variants[static restrict 1]);


bool
gva_variant_eq(GVA_Variant const lhs, GVA_Variant const rhs);


GVA_Variant
gva_variant_dup(GVA_Allocator const allocator, GVA_Variant const variant);


size_t
gva_variant_length(GVA_Variant const variant);


GVA_String
gva_patch(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const n, GVA_Variant const variants[static restrict n]);


#endif // GVA_VARIANT_H
