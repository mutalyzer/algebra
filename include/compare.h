#ifndef GVA_COMPARE_H
#define GVA_COMPARE_H


#include <stddef.h>     // size_t

#include "allocator.h"  // GVA_Allocator
#include "lcs_graph.h"  // GVA_LCS_Graph
#include "relations.h"  // GVA_Relation
#include "types.h"      // gva_uint
#include "variant.h"    // GVA_VARIANT


GVA_Relation
gva_compare_supremals(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref],
    GVA_Variant const lhs, GVA_Variant const rhs);


// FIXME: move to lcs_graph as `uniq_atomics`.
void
bitset_fill(GVA_LCS_Graph const graph,
    gva_uint const offset,
    gva_uint const start, gva_uint const end,
    size_t dels[static restrict 1],
    size_t as[static restrict 1],
    size_t cs[static restrict 1],
    size_t gs[static restrict 1],
    size_t ts[static restrict 1]);


#endif // GVA_COMPARE_H
