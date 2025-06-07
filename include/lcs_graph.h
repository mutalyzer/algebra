#ifndef GVA_LCS_GRAPH_H
#define GVA_LCS_GRAPH_H


#include <stdbool.h>    // bool
#include <stddef.h>     // size_t

#include "allocator.h"  // GVA_Allocator
#include "types.h"      // gva_uint


// Internal: Each edge tracks its tail as an index to a node (in the
// nodes array in a graph) as well as a possible `next` edge; a
// singly linked list for all edges that belong to a particular node.
typedef struct
{
    gva_uint tail;
    gva_uint next;
} GVA_Edge;


// Internal: The triple (`row`, `col`, `length`) uniquely defines each
// node in the graph. Its outgoing edges are found in the edges array of
// graph as a singly linked list with `edges` as entry point.
// A possible `lambda` edge is given as an index to another node.
typedef struct
{
    gva_uint row;
    gva_uint col;
    gva_uint length;
    gva_uint edges;
    gva_uint lambda;
} GVA_Node;


// The LCS graph stores arrays of nodes and edges with `source` as the
// entry point in the `nodes` array.
typedef struct
{
    GVA_Node* nodes;
    GVA_Edge* edges;
    GVA_Node* local_supremal;
    gva_uint  source;
} GVA_LCS_Graph;


typedef struct
{
    gva_uint start;
    gva_uint end;
    gva_uint obs_start;
    gva_uint obs_end;
} GVA_Variant;


GVA_LCS_Graph
gva_lcs_graph_init(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs],
    size_t const shift);


void
gva_lcs_graph_destroy(GVA_Allocator const allocator, GVA_LCS_Graph self);


gva_uint
gva_edges(GVA_Node const head, GVA_Node const tail,
    bool const is_source, bool const is_sink,
    GVA_Variant* const variant);


#endif // GVA_LCS_GRAPH_H
