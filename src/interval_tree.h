#ifndef GVA_INTERVAL_TREE_H
#define GVA_INTERVAL_TREE_H


#include "../include/allocator.h"   // GVA_Allocator
#include "../include/types.h"       // gva_uint


typedef struct
{
    gva_uint child[2];
    gva_uint start;
    gva_uint end;
    gva_uint max;
    // FIXME: we assume gva_uint is 32bits
    gva_uint balance   :  3;
    gva_uint inserted  : 29;
    gva_uint allele;
    gva_uint distance;
} Interval_Tree_Node;


typedef struct
{
    Interval_Tree_Node* nodes;
    gva_uint            root;
} Interval_Tree;


Interval_Tree
interval_tree_init(void);


void
interval_tree_destroy(GVA_Allocator const allocator, Interval_Tree self[static 1]);


gva_uint
interval_tree_insert(Interval_Tree self[static 1], gva_uint const idx);


#endif  // GVA_INTERVAL_TREE_H
