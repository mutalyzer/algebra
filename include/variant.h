#ifndef GVA_VARIANT_H
#define GVA_VARIANT_H


#include <stddef.h>     // size_t


#include "types.h"      // GVA_String


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


#endif // GVA_VARIANT_H
