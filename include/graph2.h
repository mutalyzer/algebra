#include <stdint.h>     // uint32_t
#include "../include/edit.h"            // VA_LCS_Node
#include "../include/variant.h"         // VA_Variant


typedef struct
{
    uint32_t tail;
    uint32_t next;
} Edge2;


typedef struct
{
    uint32_t row;
    uint32_t col;
    uint32_t length;
    uint32_t edges;
    uint32_t lambda;
} Node2;


typedef struct
{
    Node2* nodes;
    Edge2* edges;
    uint32_t source;
} Graph2;


Graph2
build(size_t const len_ref, char const reference[static len_ref],
      size_t const len_obs, char const observed[static len_obs],
      size_t const shift);
