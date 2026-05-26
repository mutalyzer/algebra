#ifndef GVA_TYPES_H
#define GVA_TYPES_H


#include <inttypes.h>   // PRIu32
#include <stdint.h>     // uint32_t, UINT32_MAX


// Print format specifier for the `gva_uint` type.
#define GVA_UINT_FMT "%" PRIu32


// An unsigned integer type to store any position on a reference
// sequence. For human DNA a 32-bit value is sufficient.
// Smaller integer types reduce storage requirements.
typedef uint32_t gva_uint;


// The maximum value for the `gva_uint` type. Used as NIL pointer in
// indexed structures as well as for unused/invalid elements.
static gva_uint const GVA_NULL = UINT32_MAX;


typedef struct
{
    gva_uint row;
    gva_uint col;
    gva_uint length;
} GVA_Match;


#endif // GVA_TYPES_H
