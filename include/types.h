#ifndef GVA_TYPES_H
#define GVA_TYPES_H


#include <stdint.h>     // uint32_t, UINT32_MAX


// An unsigned integer type to store any position on a reference
// sequence. For human DNA a 32-bit value is enough.
// Smaller integer types reduce storage requirements.
typedef uint32_t gva_uint;


// The maximum value for the `gva_uint` type.
static gva_uint const GVA_NULL = UINT32_MAX;


#endif // GVA_TYPES_H
