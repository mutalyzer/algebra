#ifndef GVA_EDIT_H
#define GVA_EDIT_H


#include <stddef.h>     // size_t

#include "allocator.h"  // GVA_Allocator
#include "types.h"      // gva_uint


typedef struct
{
    struct Match
    {
        gva_uint row;
        gva_uint col;
    }*     matches;
    char*  uniq;
    size_t max_lcs_pos;
    size_t distance;
} GVA_Matches;


size_t
gva_edit_distance(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs]);


GVA_Matches
gva_edit_distance_matches(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs]);


#endif // GVA_EDIT_H
