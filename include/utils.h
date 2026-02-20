// NOT FREESTANDING
#ifndef GVA_UTILS_H
#define GVA_UTILS_H


#include <stddef.h>     // size_t
#include <stdio.h>      // FILE

#include "allocator.h"  // GVA_Allocator
#include "lcs_graph.h"  // GVA_LCS_Graph
#include "string.h"     // GVA_String


// Reads the first sequence from a stream in FASTA format. A single header
// line (starting with '>') is expected after which the sequence is taken
// as is only omitting '\n' characters until EOF or the next '>' at the
// beginning of a line.
// The result must be deallocated by the caller (e.g. `gva_string_destroy`).
GVA_String
gva_fasta_sequence(GVA_Allocator const allocator, FILE* const stream);


void
gva_lcs_graph_dot(FILE* const stream, GVA_LCS_Graph const graph);


#endif // GVA_UTILS_H
