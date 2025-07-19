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
    // TODO: inserted seq, data
} GVA_Stabbing_Entry;


typedef struct
{
    size_t              len_ref;
    gva_uint*           start_table;
    GVA_Stabbing_Entry* entries;
} GVA_Stabbing_Index;


GVA_Stabbing_Index
gva_stabbing_index_init(GVA_Allocator const allocator, size_t const len_ref);


size_t
gva_stabbing_index_add(GVA_Allocator const allocator, GVA_Stabbing_Index index[static 1],
    size_t const start, size_t const end);


void
gva_stabbing_index_build(GVA_Allocator const allocator, GVA_Stabbing_Index const index);


void
gva_stabbing_index_destroy(GVA_Allocator const allocator, GVA_Stabbing_Index index[static 1]);


gva_uint*
gva_stabbing_index_self_intersect(GVA_Allocator const allocator, GVA_Stabbing_Index const index);


gva_uint*
gva_stabbing_index_intersect(GVA_Allocator const allocator, GVA_Stabbing_Index const index,
    gva_uint const start, gva_uint const end);


#endif  // GVA_STABBING_H
