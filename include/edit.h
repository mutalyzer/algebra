#ifndef GVA_EDIT_H
#define GVA_EDIT_H


#include <stddef.h>     // size_t

#include "allocator.h"  // GVA_Allocator


size_t
gva_edit_distance(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs]);


#endif // GVA_EDIT_H
