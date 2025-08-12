#ifndef GVA_STABBING_H
#define GVA_STABBING_H


#include <stddef.h>     // size_t

#include "allocator.h"      // GVA_Allocator
#include "types.h"          // gva_uint


typedef struct
{
    gva_uint parent;
    gva_uint leftsibling;
    gva_uint rightchild;
    gva_uint start;
    gva_uint end;
    gva_uint inserted;
    gva_uint distance;
    gva_uint allele;
    gva_uint next;
} GVA_Stabbing_Entry;


typedef struct
{
    size_t   data;
    gva_uint distance;
    gva_uint head;
} GVA_Stabbing_Allele;


// FIXME: include trie as part of the index.
typedef struct
{
    size_t               len_ref;
    gva_uint*            start_table;
    GVA_Stabbing_Entry*  entries;
    GVA_Stabbing_Allele* alleles;
} GVA_Stabbing_Index;


GVA_Stabbing_Index
gva_stabbing_index_init(GVA_Allocator const allocator, size_t const len_ref);


void
gva_stabbing_index_destroy(GVA_Allocator const allocator, GVA_Stabbing_Index index[static 1]);


size_t
gva_stabbing_index_add_allele(GVA_Allocator const allocator, GVA_Stabbing_Index index[static 1],
    size_t const data, gva_uint const distance);


size_t
gva_stabbing_index_add_part(GVA_Allocator const allocator, GVA_Stabbing_Index index[static 1],
    gva_uint const start, gva_uint const end, gva_uint const inserted,
    gva_uint const distance, gva_uint const allele);


void
gva_stabbing_index_build(GVA_Allocator const allocator, GVA_Stabbing_Index const index);


gva_uint*
gva_stabbing_index_self_intersect(GVA_Allocator const allocator, GVA_Stabbing_Index const index);


gva_uint*
gva_stabbing_index_intersect(GVA_Allocator const allocator, GVA_Stabbing_Index const index,
    gva_uint const start, gva_uint const end);


#endif  // GVA_STABBING_H
