#ifndef GVA_LCS_GRAPH_H
#define GVA_LCS_GRAPH_H


#include <stdbool.h>    // bool
#include <stddef.h>     // size_t

#include "allocator.h"  // GVA_Allocator
#include "string.h"     // GVA_String
#include "types.h"      // GVA_Match, gva_uint
#include "variant.h"    // GVA_Variant


// Internal: Each edge tracks its tail as an index to a node (in the
// nodes array in a graph) as well as a possible `next` edge; a
// singly linked list for all edges that belong to a particular node.
typedef struct
{
    gva_uint tail;
    gva_uint next;
} GVA_Edge;


// Internal: The `GVA_Match` structure uniquely defines each
// node in the graph. Its outgoing edges are found in the edges array of
// graph as a singly linked list with `edges` as entry point.
// A possible `lambda` edge is given as an index to another node.
typedef struct
{
    GVA_Match match;
    gva_uint  edges;
    gva_uint  lambda;
} GVA_Node;


typedef struct
{
    GVA_Match match;
    gva_uint  distance;
    gva_uint  link;
} GVA_Dom_Node;


// The LCS graph stores arrays of nodes and edges with `source` as the
// entry point in the `nodes` array.
typedef struct
{
    GVA_Node*     nodes;
    GVA_Edge*     edges;
    GVA_Dom_Node* dom_nodes;
    GVA_String    observed;
    gva_uint      source;
} GVA_LCS_Graph;


GVA_LCS_Graph
gva_lcs_graph_init(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const len_obs, char const observed[static restrict len_obs],
    size_t const offset);


void
gva_lcs_graph_destroy(GVA_Allocator const allocator,
    GVA_LCS_Graph self, bool const observed);


GVA_LCS_Graph
gva_lcs_graph_from_allele(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const n, GVA_Variant const variants[static restrict n]);


GVA_LCS_Graph
gva_lcs_graph_from_variants(GVA_Allocator const allocator,
    size_t const len_ref, char const reference[static restrict len_ref],
    size_t const n, GVA_Variant const variants[static restrict n]);


void
gva_lcs_graph_uniq_atomics(GVA_LCS_Graph const self,
    size_t const offset,
    size_t const start, size_t const end,
    size_t dels[static restrict 1],
    size_t as[static restrict 1],
    size_t cs[static restrict 1],
    size_t gs[static restrict 1],
    size_t ts[static restrict 1]);


size_t
gva_lcs_graph_distance(GVA_LCS_Graph const self);


GVA_Variant
gva_lcs_graph_local_supremal(GVA_LCS_Graph const self,
    size_t const start, size_t const end);


GVA_Variant
gva_lcs_graph_supremal(GVA_LCS_Graph const self);


size_t
gva_edges(char const observed[static restrict 1],
    GVA_Match const head, GVA_Match const tail,
    bool const is_source, bool const is_sink,
    GVA_Variant variant[static restrict 1]);


#endif // GVA_LCS_GRAPH_H
