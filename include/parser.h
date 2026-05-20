#ifndef GVA_PARSER_H
#define GVA_PARSER_H


#include <stddef.h>     // size

#include "allocator.h"  // GVA_Allocator
#include "variant.h"    // GVA_Variant


GVA_Variant*
gva_parse_hgvs(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len, char const expression[static restrict len]);


size_t
gva_parse_spdi(size_t const len, char const expression[static restrict len],
    GVA_Variant variants[static restrict 1]);


#endif  // GVA_PARSER_H
