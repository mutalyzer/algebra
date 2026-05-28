#ifndef GVA_INTERVAL_TREE_H
#define GVA_INTERVAL_TREE_H


#include <limits.h>     // CHAR_BIT
#include <stdint.h>     // int8_t

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/types.h"       // gva_uint


// [-4, ..., 3], we use [-2, ..., 2]
#define BALANCE_BITS 3


typedef struct
{
    gva_uint child[2];
    gva_uint start;
    gva_uint end;
    gva_uint max;
    int8_t   balance  : BALANCE_BITS;
    gva_uint inserted : sizeof(gva_uint) * CHAR_BIT - BALANCE_BITS;
    gva_uint alleles;
    gva_uint distance;
    gva_uint dfa;
} Interval_Tree_Node;


typedef struct
{
    Interval_Tree_Node* nodes;
    GVA_Allocator       allocator;
    gva_uint            root;
} Interval_Tree;


Interval_Tree
interval_tree_init(GVA_Allocator const allocator, gva_uint const root);


void
interval_tree_destroy(Interval_Tree self[static 1]);


gva_uint
interval_tree_insert(Interval_Tree self[static 1], gva_uint const idx);


gva_uint*
interval_tree_intersection(GVA_Allocator const allocator, Interval_Tree const self,
    gva_uint const start, gva_uint const end);


#endif  // GVA_INTERVAL_TREE_H
