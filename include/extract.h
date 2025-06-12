#ifndef GVA_EXTRACT_H
#define GVA_EXTRACT_H


#include "allocator.h"      // GVA_Allocator
#include "lcs_graph.h"      // GVA_LCS_Graph, GVA_Variant


GVA_Variant*
gva_extract(GVA_Allocator const allocator, GVA_LCS_Graph const graph);


#endif // GVA_EXTRACT_H
