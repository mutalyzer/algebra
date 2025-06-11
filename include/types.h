#ifndef GVA_TYPES_H
#define GVA_TYPES_H


#include <stdint.h>     // uint32_t, UINT32_MAX


#define GVA_VARIANT_PRINT(variant) variant.start, variant.end, (int) variant.sequence.len, variant.sequence.str
#define GVA_VARIANT_FMT "%u:%u/%.*s"


// An unsigned integer type to store any position on a reference
// sequence. For human DNA a 32-bit value is enough.
// Smaller integer types reduce storage requirements.
typedef uint32_t gva_uint;


// The maximum value for the `gva_uint` type.
static gva_uint const GVA_NULL = UINT32_MAX;


typedef struct
{
    char const* str;
    gva_uint    len;
} GVA_String;


#endif // GVA_TYPES_H
