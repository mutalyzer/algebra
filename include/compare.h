#ifndef GVA_COMPARE_H
#define GVA_COMPARE_H


#include <stddef.h>     // size_t

#include "allocator.h"  // GVA_Allocator
#include "lcs_graph.h"  // GVA_LCS_Graph
#include "relations.h"  // GVA_Relation


GVA_Relation
gva_compare_graphs(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref],
    GVA_LCS_Graph const lhs, GVA_LCS_Graph const rhs);


#endif // GVA_COMPARE_H
