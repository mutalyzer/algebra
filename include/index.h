#ifndef GVA_INDEX_H
#define GVA_INDEX_H


#include <stddef.h>     // size_t

#include "allocator.h"  // GVA_Allocator
#include "types.h"      // gva_uint
#include "variant.h"    // GVA_Variant


typedef struct
{
    GVA_Allocator       allocator;
    gva_uint            root;
    struct Allele*      alleles;
    struct Node*        nodes;
    struct Next_Allele* next_alleles;
} GVA_Index;


GVA_Index
gva_index_init(GVA_Allocator const allocator);


void
gva_index_destroy(GVA_Index self);


void
gva_index_add(GVA_Index self, size_t const data,
    size_t const n, GVA_Variant const variants[static n]);


#endif  // GVA_INDEX_H

