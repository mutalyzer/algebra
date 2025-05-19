#ifndef EDIT_H
#define EDIT_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>     // size_t
#include <stdint.h>     // uint32_t

#include "alloc.h"      // VA_Allocator


typedef struct
{
    uint32_t row;
    uint32_t col;
    uint32_t length;
    uint32_t incoming;
    uint32_t idx;
    uint32_t x;
} VA_LCS_Node;


typedef struct
{
    uint32_t row;
    uint32_t col;
    uint32_t length;
    uint32_t lcs_pos;
    uint32_t prev;
    uint32_t idx;
} VA_LCS_Node2;


size_t
va_edit_distance_only(VA_Allocator const allocator,
                      size_t const len_ref,
                      char const reference[static restrict len_ref],
                      size_t const len_obs,
                      char const observed[static restrict len_obs]);


size_t
va_edit(VA_Allocator const allocator,
        size_t const len_ref,
        char const reference[static restrict len_ref],
        size_t const len_obs,
        char const observed[static restrict len_obs],
        VA_LCS_Node*** restrict lcs_nodes);


size_t
va_edit2(VA_Allocator const allocator,
         size_t const len_ref,
         char const reference[static restrict len_ref],
         size_t const len_obs,
         char const observed[static restrict len_obs],
         VA_LCS_Node2** restrict lcs_nodes,
         uint32_t** restrict lcs_index);


#ifdef __cplusplus
} // extern "C"
#endif

#endif
