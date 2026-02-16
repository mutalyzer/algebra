#ifndef GVA_INDEX_H
#define GVA_INDEX_H


#include <stddef.h>     // size_t

#include "allocator.h"  // GVA_Allocator
#include "lcs_graph.h"  // GVA_LCS_Graph
#include "relations.h"  // GVA_Relation
#include "variant.h"    // GVA_Variant


typedef struct GVA_Index GVA_Index;  // opaque


typedef struct
{
    GVA_String   allele;
    GVA_Relation relation;
} GVA_Result;


GVA_Index*
gva_index_init(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static len_ref]);


GVA_Index*
gva_index_destroy(GVA_Index* const self);


// insertions need to happen in allele-order
void
gva_index_insert(GVA_Index* restrict const self,
    size_t const len_id, char const id[static restrict len_id],
    GVA_Variant const variant, size_t const distance);


// FIXME: returns
GVA_Result*
gva_index_query(GVA_Allocator const allocator,
    GVA_Index* const self, GVA_LCS_Graph const graph);


#endif // GVA_INDEX_H
