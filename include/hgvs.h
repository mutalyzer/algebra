#ifndef GVA_HGVS_H
#define GVA_HGVS_H


#include <stddef.h>     // size_t

#include "variant.h"    // GVA_Variant


size_t
gva_to_hgvs(size_t const len, char buffer[static restrict len],
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const n, GVA_Variant const variants[static restrict n]);


#endif // GVA_HGVS_H
