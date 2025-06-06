#ifndef GVA_ALIGN_H
#define GVA_ALIGN_H


#include <stddef.h>     // size_t

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/types.h"       // gva_uint


typedef struct
{
    gva_uint row;
    gva_uint col;
    gva_uint length : 31;
    gva_uint moved  :  1;
    gva_uint incoming;
    gva_uint idx;
    gva_uint next;
} LCS_Node;


typedef struct
{
    size_t length;
    struct
    {
        gva_uint head;
        gva_uint tail;
    }* index;
    LCS_Node* nodes;
} LCS_Alignment;


LCS_Alignment
lcs_align(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs]);


#endif // GVA_ALIGN_H
