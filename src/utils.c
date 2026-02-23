// NOT FREESTANDING
#include <errno.h>      // errno
#include <stddef.h>     // NULL, size_t
#include <stdlib.h>     // rand
#include <stdio.h>      // FILE, stderr, fgets, fprintf
#include <string.h>     // strcspn

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_edges
#include "../include/string.h"      // GVA_String
#include "../include/utils.h"       // gva_fasta_sequence
#include "array.h"      // array_length


GVA_String
gva_fasta_sequence(GVA_Allocator const allocator, FILE* const stream)
{
    static size_t const FASTA_LINE_SIZE = 65536;
    size_t capacity = FASTA_LINE_SIZE;
    GVA_String sequence = {0, allocator.allocate(allocator.context, NULL, 0, capacity)};
    if (sequence.str == NULL)
    {
        return (GVA_String) {0};
    } // if

    while (fgets((char*) sequence.str + sequence.len, FASTA_LINE_SIZE, stream) != NULL)
    {
        if (sequence.str[sequence.len] == '>')
        {
            if (sequence.len == 0)
            {
                continue;
            } // if
            break;
        } // if
        sequence.len += strcspn(sequence.str + sequence.len, "\n");

        if (capacity - FASTA_LINE_SIZE < sequence.len)
        {
            sequence.str = allocator.allocate(allocator.context, (char*) sequence.str, capacity, capacity + FASTA_LINE_SIZE);
            if (sequence.str == NULL)
            {
                sequence.len = 0;
                return sequence;
            } // if
            capacity += FASTA_LINE_SIZE;
        } // if
    } // while
    sequence.str = allocator.allocate(allocator.context, (char*) sequence.str, capacity, sequence.len);
    return sequence;
} // gva_fasta_sequence


GVA_String
gva_fasta_sequence_blob(GVA_Allocator const allocator, FILE* const stream)
{
    GVA_String reference = {0};

    errno = 0;
    if (fread(&reference.len, sizeof(reference.len), 1, stream) != 1)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return (GVA_String) {0};
    }  // if

    reference.str = allocator.allocate(allocator.context, NULL, 0, reference.len);
    if (reference.str == NULL)
    {
        fprintf(stderr, "OOM\n");
        return (GVA_String) {0};
    } // if

    errno = 0;
    if (fread((char*) reference.str, 1, reference.len, stream) != reference.len)
    {
        fprintf(stderr, "error: %s\n", strerror(errno));
        gva_string_destroy(allocator, reference);
        return (GVA_String) {0};
    } // if

    return reference;
} // gva_fasta_sequence_blob


void
gva_lcs_graph_dot(FILE* const stream, GVA_LCS_Graph const graph)
{
    fprintf(stream, "strict digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\ni[label=\"\",shape=none,width=0]\ni->%u\n", graph.source);
    for (size_t i = 0; i < array_length(graph.nodes); ++i)
    {
        fprintf(stream, "%zu[label=\"(%u, %u, %u)\"%s]\n", i, graph.nodes[i].match.row, graph.nodes[i].match.col, graph.nodes[i].match.length, graph.nodes[i].edges == GVA_NULL ? ",peripheries=2" : "");
        if (graph.nodes[i].lambda != GVA_NULL)
        {
            fprintf(stream, "%zu->%u[label=\"&lambda;\",style=dashed]\n", i, graph.nodes[i].lambda);
        } // if
        for (gva_uint j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
        {
            GVA_Variant variant;
            gva_uint const count = gva_edges(graph.observed.str,
                graph.nodes[i].match, graph.nodes[graph.edges[j].tail].match,
                i == graph.source, graph.nodes[graph.edges[j].tail].edges == GVA_NULL,
                &variant);
            if (count > 1)
            {
                fprintf(stream, "%zu->%u[label=\"" GVA_VARIANT_FMT " x %u\",penwidth=2]\n", i, graph.edges[j].tail, GVA_VARIANT_PRINT(variant), count);
            } // if
            else
            {
                fprintf(stream, "%zu->%u[label=\"" GVA_VARIANT_FMT "\"]\n", i, graph.edges[j].tail, GVA_VARIANT_PRINT(variant));
            } // else
        } // for
    } // for
    for (size_t i = 0; i < array_length(graph.dom_nodes); ++i)
    {
        fprintf(stream, "%u[penwidth=2]\n", graph.dom_nodes[i].link);
    } // for
    fprintf(stream, "}\n");
} // gva_lcs_graph_dot
