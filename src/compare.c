#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint8_t

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/compare.h"     // gva_compare*
#include "../include/edit.h"        // gva_edit_distance
#include "../include/relations.h"   // GVA_Relation, GVA_CONTAINS, GVA_DISJOINT,
                                    // GVA_EQUIVALENT GVA_IS_CONTAINED, GVA_OVERLAP
#include "../include/variant.h"     // GVA_Variant, gva_patch, gva_variant_eq
#include "common.h"     // MAX, MIN
#include "dfa.h"        // dfa_overlap


size_t
gva_compare_distance(GVA_Allocator const allocator,
    size_t const len, char const reference[static len],
    GVA_Variant const lhs, GVA_Variant const rhs)
{
    size_t const start = MIN(lhs.start, rhs.start);
    size_t const end = MAX(lhs.end, rhs.end);

    size_t const len_lhs = (lhs.start - start) + lhs.sequence.len + (end - lhs.end);
    size_t const len_rhs = (rhs.start - start) + rhs.sequence.len + (end - rhs.end);

    if (len_lhs == 0)
    {
        return len_rhs;
    } // if
    if (len_rhs == 0)
    {
        return len_lhs;
    } // if

    GVA_String observed_lhs = gva_patch(allocator, end - start, reference + start,
        1, &(GVA_Variant const) {lhs.start - start, lhs.end - start, lhs.sequence});
    GVA_String observed_rhs = gva_patch(allocator, end - start, reference + start,
        1, &(GVA_Variant const) {rhs.start - start, rhs.end - start, rhs.sequence});

    size_t const distance = gva_edit_distance(allocator,
        observed_lhs.len, observed_lhs.str, observed_rhs.len, observed_rhs.str);

    gva_string_destroy(allocator, observed_rhs);
    gva_string_destroy(allocator, observed_lhs);

    return distance;
} // gva_compare_distance


size_t
gva_compare_included(GVA_Allocator const allocator,
    size_t const len, char const reference[static restrict len],
    GVA_Variant const lhs, size_t const lhs_distance, uint8_t lhs_dfa[static restrict 1],
    GVA_Variant const rhs, size_t const rhs_distance, uint8_t rhs_dfa[static restrict 1])
{
    size_t const distance = gva_compare_distance(allocator, len, reference, lhs, rhs);

    if (distance == 0 || rhs_distance - lhs_distance == distance)
    {
        return lhs_distance;
    } // if
    if (lhs_distance - rhs_distance == distance)
    {
        return rhs_distance;
    } // if

    return dfa_overlap(allocator, len, reference, lhs, lhs_dfa, rhs, rhs_dfa);
} // gva_compare_included


inline GVA_Relation
gva_compare_relation_from_included(size_t const included,
    size_t const lhs_distance, size_t const rhs_distance,
    size_t excluded[static 1])
{
    *excluded = lhs_distance + rhs_distance - 2 * included;
    if (*excluded == 0)
    {
        return GVA_EQUIVALENT;
    } // if
    if (included == 0)
    {
        return GVA_DISJOINT;
    } // if
    if (lhs_distance == included)
    {
        return GVA_IS_CONTAINED;
    } // if
    if (rhs_distance == included)
    {
         return GVA_CONTAINS;
    } // if
    return GVA_OVERLAP;
} // gva_compare_relation_from_included


inline GVA_Relation
gva_compare(GVA_Allocator const allocator,
    size_t const len, char const reference[static restrict len],
    GVA_Variant const lhs, size_t const lhs_distance, uint8_t lhs_dfa[static restrict 1],
    GVA_Variant const rhs, size_t const rhs_distance, uint8_t rhs_dfa[static restrict 1])
{
    if (gva_variant_eq(lhs, rhs))
    {
        return GVA_EQUIVALENT;
    } // if
    size_t const included = gva_compare_included(allocator, len, reference,
        lhs, lhs_distance, lhs_dfa,
        rhs, rhs_distance, rhs_dfa);
    return gva_compare_relation_from_included(included, lhs_distance, rhs_distance, &(size_t) {0});
} // gva_compare
