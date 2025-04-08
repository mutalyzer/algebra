#include "../include/alloc.h"           // VA_Allocator
#include "../include/graph.h"           // Graph


void
local_supremal(VA_Allocator const allocator, Graph const graph, size_t const len_obs, char const observed[static len_obs], bool const debug);

void
canonical(VA_Allocator const allocator, Graph const graph, size_t const len_obs, char const observed[static len_obs], bool const debug);
