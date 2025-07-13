#ifndef GVA_COMPARE_H
#define GVA_COMPARE_H


#include "allocator.h"  // GVA_Allocator
#include "relations.h"  // GVA_Relation
#include "variant.h"    // GVA_VARIANT


extern size_t disjoint_count;


GVA_Relation
gva_compare(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref],
    GVA_Variant const lhs, GVA_Variant const rhs);


#endif // GVA_COMPARE_H
