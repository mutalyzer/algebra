#ifndef GVA_EXTRACT_H
#define GVA_EXTRACT_H


#include <stddef.h>     // size_t


#include "allocator.h"      // GVA_Allocator
#include "lcs_graph.h"      // GVA_LCS_Graph


void
gva_extract(GVA_Allocator const allocator,
    GVA_LCS_Graph const graph, size_t const len_obs, char const observed[static len_obs]);


#endif // GVA_EXTRACT_H
