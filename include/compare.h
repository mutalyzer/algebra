#ifndef GVA_COMPARE_H
#define GVA_COMPARE_H


#include <stddef.h>     // size_t
#include <stdint.h>     // uint8_t

#include "allocator.h"  // GVA_Allocator
#include "relations.h"  // GVA_Relation
#include "variant.h"    // GVA_Variant


size_t
gva_compare_distance(GVA_Allocator const allocator,
    size_t const len, char const reference[static len],
    GVA_Variant const lhs, GVA_Variant const rhs);


size_t
gva_compare_included(GVA_Allocator const allocator,
    size_t const len, char const reference[static restrict len],
    GVA_Variant const lhs, size_t const lhs_distance, uint8_t lhs_dfa[static restrict 1],
    GVA_Variant const rhs, size_t const rhs_distance, uint8_t rhs_dfa[static restrict 1]);


GVA_Relation
gva_compare_relation_from_included(size_t const included,
    size_t const lhs_distance, size_t const rhs_distance,
    size_t excluded[static 1]);


GVA_Relation
gva_compare(GVA_Allocator const allocator,
    size_t const len, char const reference[static restrict len],
    GVA_Variant const lhs, size_t const lhs_distance, uint8_t lhs_dfa[static restrict 1],
    GVA_Variant const rhs, size_t const rhs_distance, uint8_t rhs_dfa[static restrict 1]);


#endif // GVA_COMPARE_H
