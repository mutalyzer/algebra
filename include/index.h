#ifndef GVA_INDEX_H
#define GVA_INDEX_H


#include <stddef.h>     // size_t

#include "allocator.h"  // GVA_Allocator
#include "lcs_graph.h"  // GVA_LCS_Graph
#include "relations.h"  // GVA_Relation
#include "types.h"      // GVA_Interval, gva_uint
#include "variant.h"    // GVA_Variant


typedef struct GVA_Index GVA_Index;  // opaque


typedef struct
{
    struct GVA_Query_Allele
    {
        gva_uint     idx;
        GVA_Relation relation;
        gva_uint     included;
        gva_uint     excluded;
        GVA_Interval hits;
    }* alleles;
    struct GVA_Query_Hit
    {
        GVA_Relation relation;
        gva_uint     included;
        gva_uint     excluded;
        GVA_Interval index;
        GVA_Interval query;
    }* hits;
} GVA_Query_Result;


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


GVA_Query_Result
gva_index_query(GVA_Allocator const allocator,
    GVA_Index* const self, GVA_LCS_Graph const graph);


GVA_String
gva_index_id(GVA_Index const* const self, size_t const idx);


GVA_Variant
gva_index_variant(GVA_Index const* const self,
    size_t const allele_idx, size_t const variant_idx);


#endif // GVA_INDEX_H
