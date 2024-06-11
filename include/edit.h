#ifndef EDIT_H
#define EDIT_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>     // size_t

#include "alloc.h"      // VA_Allocator


typedef struct
{
    size_t row;
    size_t col;
    size_t length;
} VA_LCS_Node;


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


#ifdef __cplusplus
} // extern "C"
#endif

#endif
