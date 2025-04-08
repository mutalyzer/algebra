#include <stdint.h>     // uint32_t
#include "../include/edit.h"            // VA_LCS_Node
#include "../include/variant.h"         // VA_Variant

typedef struct
{
    uint32_t tail;
    uint32_t next;
    VA_Variant variant;
} Edge;


typedef struct
{
    uint32_t row;
    uint32_t col;
    uint32_t length;
    uint32_t edges;
    uint32_t lambda;
} Node;

typedef struct
{
    Node* nodes;
    Edge* edges;
    VA_Variant supremal;
    uint32_t source;
} Graph;


Graph
build_graph(VA_Allocator const allocator,
            size_t const len_ref, size_t const len_obs,
            size_t const len_lcs, VA_LCS_Node* lcs_nodes[static len_lcs],
            size_t const shift, bool const debug);


void
to_dot(Graph const graph, size_t const len_obs, char const observed[static len_obs]);


void
to_json(Graph const graph, size_t const len_obs, char const observed[static len_obs], bool const lambda);


size_t
edges(uint32_t const head_row, uint32_t const head_col, uint32_t const head_length,
      uint32_t const tail_row, uint32_t const tail_col, uint32_t const tail_length,
      bool const is_source,
      size_t const len_obs, char const observed[static len_obs],
      VA_Variant const edge);


void
destroy(VA_Allocator const allocator, Graph* const graph);
