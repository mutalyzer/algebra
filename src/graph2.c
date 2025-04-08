#include <stdbool.h>    // bool, false, true
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // intmax_t, uint32_t
#include <stdio.h>      // stderr, fprintf, printf

#include "../include/alloc.h"           // VA_Allocator
#include "../include/array.h"           // va_array_*
#include "../include/edit.h"            // VA_LCS_Node, va_edit
#include "../include/std_alloc.h"       // va_std_allocator
#include "../include/variant.h"         // VA_Variant
#include "../include/graph2.h"          // Edge2, Graph2, Node2


#define MAX(lhs, rhs) (((lhs) > (rhs)) ? (lhs) : (rhs))
#define MIN(lhs, rhs) (((lhs) < (rhs)) ? (lhs) : (rhs))

#define print_variant(variant, observed) variant.start, variant.end, (int) variant.obs_end - variant.obs_start, observed + variant.obs_start
#define VAR_FMT "%u:%u/%.*s"


static uint32_t const GVA_NULL = UINT32_MAX;



static uint32_t
edges2(VA_LCS_Node const head, VA_LCS_Node const tail,
       bool const is_source, bool const is_sink,
       size_t const len_obs, char const observed[static len_obs], VA_Variant* const variant)
{
//    printf("%u %u %u -> %u %u %u %s %s\n", head.row, head.col, head.length, tail.row, tail.col, tail.length, is_source ? "SOURCE" : "", is_sink ? "SINK" : "");

    intmax_t const row = (intmax_t) head.row - is_source;
    intmax_t const col = (intmax_t) head.col - is_source;
    uint32_t const head_length = head.length + is_source;
    uint32_t const tail_length = tail.length + is_sink;

    intmax_t const offset = MIN((intmax_t) tail.row - row, (intmax_t) tail.col - col) - 1;

    uint32_t const head_offset = offset > 0 ? MIN(head_length, offset + 1) : 1;
    uint32_t const tail_offset = offset < 0 ? MIN(tail_length, -offset) : 0;

    //printf("    %ld: (%u, %u)\n", offset, head_offset, tail_offset);

    if (head_offset > head_length || (tail_length > 0 && tail_offset >= tail_length))
    {
        fprintf(stderr, "NO EDGE\n");
        return 0;
    } // if

    uint32_t const count = MIN(head_length - head_offset, tail_length - tail_offset - 1) + 1;

    *variant = (VA_Variant) {row + head_offset, tail.row + tail_offset, col + head_offset, tail.col + tail_offset};

    fprintf(stderr, VAR_FMT " x %u\n", print_variant((*variant), observed), count);
    return count;
} // edges2


static size_t
bfs_traversal2(Graph2 const graph, size_t const len_obs, char const observed[static len_obs])
{
    struct
    {
        uint32_t depth;
        uint32_t next;
    }* table = va_std_allocator.alloc(va_std_allocator.context, NULL, 0, sizeof(*table) * va_array_length(graph.nodes));

    if (table == NULL)
    {
        return -1;
    } // if

    for (size_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        table[i].depth = 0;
        table[i].next = GVA_NULL;
    } // for
    uint32_t head = graph.source;
    uint32_t tail = graph.source;

    size_t count = 0;
    fprintf(stderr, "digraph{\nrankdir=LR\nedge[fontname=monospace]\nnode[fixedsize=true,fontname=serif,shape=circle,width=1]\nsi[label=\"\",shape=none,width=0]\n");
    fprintf(stderr, "si->s%u\n", graph.source);

    bool first = true;

    printf("    \"edges\": [\n");
    while (head != GVA_NULL)
    {
        uint32_t len = graph.nodes[head].length;
        //printf("pop %u\n", head);
        for (uint32_t i = head; i != GVA_NULL; i = graph.nodes[i].lambda)
        {
            fprintf(stderr, "s%u[label=\"(%u, %u, %u)\"%s]\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, graph.nodes[i].edges == GVA_NULL && graph.nodes[i].lambda == GVA_NULL ? ",peripheries=2" : "");
            if (i != head)
            {
                len += graph.nodes[i].length;
                //printf("lambda %u\n", i);
            } // if
            for (uint32_t j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
            {
                if (count > 0)
                {
                    //printf(",\n");
                } // if
                count += 1;
                VA_Variant variant;
                uint32_t const count2 = edges2(((VA_LCS_Node) {graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, 0, 0}),
                                               ((VA_LCS_Node) {graph.nodes[graph.edges[j].tail].row, graph.nodes[graph.edges[j].tail].col, graph.nodes[graph.edges[j].tail].length, 0, 0}),
                                               graph.nodes[head].row == graph.nodes[graph.source].row && graph.nodes[head].col == graph.nodes[graph.source].col, graph.nodes[graph.edges[j].tail].edges == GVA_NULL, len_obs, observed, &variant);

                for (uint32_t k = 0; k < count2; ++k)
                {
                    if (!first)
                    {
                        printf(",\n");
                    } // if
                    first = false;
                    fprintf(stderr, "s%u->s%u[label=\"" VAR_FMT "\"]\n", head, graph.edges[j].tail, print_variant(variant, observed));
                    printf("         {\"head\": \"s%u\", \"tail\": \"s%u\", \"variant\": \"" VAR_FMT "\"}", head, graph.edges[j].tail, print_variant(variant, observed));
                    variant.start += 1;
                    variant.end += 1;
                    variant.obs_start += 1;
                    variant.obs_end += 1;
                } // for

                if (table[graph.edges[j].tail].depth > 0)
                {
                    //printf("skip %u\n", graph.edges[j].tail);
                    continue;
                } // if

                //printf("push %u\n", graph.edges[j].tail);
                table[graph.edges[j].tail].depth = table[i].depth + 1;
                table[tail].next = graph.edges[j].tail;
                tail = graph.edges[j].tail;
            } // for
        } // for
        head = table[head].next;

        /*
        printf("  # \tdepth\tnext\n");
        for (size_t i = 0; i < va_array_length(graph.nodes); ++i)
        {
            printf("%3zu:\t%5u\t%4d\n", i, table[i].depth, table[i].next);
        } // for
        printf("head: %d\ntail: %d\n", head, tail);
        */
    } // while

    fprintf(stderr, "}\n");
    printf("\n    ]\n");

    table = va_std_allocator.alloc(va_std_allocator.context, table, sizeof(*table) * va_array_length(graph.nodes), 0);

    return count;
} // bfs_traversal2


Graph2
build(size_t const len_ref, char const reference[static len_ref],
      size_t const len_obs, char const observed[static len_obs],
      size_t const shift)
{
    //printf("build!\n");
    VA_LCS_Node** lcs_nodes = NULL;
    size_t const len_lcs = va_edit(va_std_allocator, len_ref, reference, len_obs, observed, &lcs_nodes);

    Graph2 graph = {NULL, NULL, GVA_NULL};

    if (len_lcs == 0 || lcs_nodes == NULL)
    {
        va_array_append(va_std_allocator, graph.nodes, ((Node2) {shift, shift, 0, GVA_NULL, GVA_NULL}));
        graph.source = va_array_length(graph.nodes) - 1;

        if (len_ref == 0 && len_obs == 0)
        {
            //graph.supremal = (VA_Variant) {0, 0, 0, 0};
            return graph;
        } // if
        va_array_append(va_std_allocator, graph.nodes, ((Node2) {len_ref, len_obs, 0, GVA_NULL, GVA_NULL}));
        uint32_t const sink = va_array_length(graph.nodes) - 1;
        //graph.supremal = (VA_Variant) {shift, shift + len_ref, 0, len_obs};

        va_array_append(va_std_allocator, graph.edges, ((Edge2) {sink, graph.nodes[graph.source].edges}));
        graph.nodes[graph.source].edges = va_array_length(graph.edges) - 1;
        return graph;
    } // if

    VA_LCS_Node* last = &lcs_nodes[len_lcs - 1][va_array_length(lcs_nodes[len_lcs - 1]) - 1];
    if (last->row + last->length != len_ref + shift || last->col + last->length != len_obs + shift)
    {
        va_array_append(va_std_allocator, lcs_nodes[len_lcs - 1],
                        ((VA_LCS_Node) {.row = len_ref + shift, .col = len_obs + shift, .length = 0}));
        last = &lcs_nodes[len_lcs - 1][va_array_length(lcs_nodes[len_lcs - 1]) - 1];
    } // if
    va_array_append(va_std_allocator, graph.nodes, ((Node2) {last->row, last->col, last->length, GVA_NULL, GVA_NULL}));
    last->idx = va_array_length(graph.nodes) - 1;
    fprintf(stderr, "MAKE NODE %u: (%u, %u, %u)\n", last->idx, last->row, last->col, last->length);
    last->incoming = -1;

    /*
    fprintf(stderr, "%zu\n", len_lcs);
    for (size_t i = 0; i < len_lcs; ++i)
    {
        fprintf(stderr, "%zu:  ", i);
        for (size_t j = 0; j < va_array_length(lcs_nodes[i]); ++j)
        {
            fprintf(stderr, "(%u, %u, %u) ", lcs_nodes[i][j].row, lcs_nodes[i][j].col, lcs_nodes[i][j].length);
        } // for
        fprintf(stderr, "\n");
    } // for
    */

    VA_LCS_Node source = {shift, shift, 0, -1, GVA_NULL};
    bool found_source = last->row == shift && last->col == shift;
    bool is_sink = true;
    for (size_t t_i = 0; t_i < len_lcs; ++t_i)
    {
        size_t const t_len = va_array_length(lcs_nodes[len_lcs - t_i - 1]);
        for (size_t t_j = 0; t_j < t_len; ++t_j)
        {
            VA_LCS_Node* const tail = &lcs_nodes[len_lcs - t_i - 1][t_len - t_j - 1];
            if (tail->idx == GVA_NULL)
            {
                // sink unreachable from tail, do not add to graph
                continue;
            } // if
            fprintf(stderr, "%u %u %u %d\n", tail->row, tail->col, tail->length, tail->incoming);

            uint32_t split_idx = tail->idx;
            for (size_t h_i = 0; h_i < MIN(len_lcs - t_i, tail->length + 1); ++h_i)
            {
                size_t const h_len = h_i == 0 ? t_len - t_j - 1: va_array_length(lcs_nodes[len_lcs - t_i - h_i - 1]);
                //printf("    level: %zu with %zu nodes\n", len_lcs - t_i - h_i - 1, h_len);
                for (size_t h_j = 0; h_j < h_len; ++h_j)
                {
                    VA_LCS_Node* head = &lcs_nodes[len_lcs - t_i - h_i - 1][h_j];

                    fprintf(stderr, "    %u %u %u: ", head->row, head->col, head->length);

                    if (head->row + head->length + h_i > tail->row + tail->length ||
                        head->col + head->length + h_i > tail->col + tail->length)
                    {
                        if (tail->row + tail->length > head->row + head->length + h_i ||
                            tail->col + tail->length > head->col + head->length + h_i)
                        {
                            fprintf(stderr, "--\n");  // FIXME: skip the remainder of this (head) level
                        } // if
                        else
                        {
                            fprintf(stderr, "CONVERSE\n");
                            VA_Variant variant;
                            if (head->idx != GVA_NULL && edges2(*tail, *head, tail->row == shift && tail->col == shift,
                                                                false, len_obs, observed, &variant))
                            {
                                va_array_append(va_std_allocator, graph.edges, ((Edge2) {head->idx, graph.nodes[tail->idx].edges}));
                                graph.nodes[tail->idx].edges = va_array_length(graph.edges) - 1;
                                fprintf(stderr, "MAKE EDGE\n");
                                tail->incoming = MIN(tail->incoming, variant.start);
                            } // if
                        } // else
                    } // if
                    else
                    {
                        VA_Variant variant;
                        uint32_t const count = edges2(*head, *tail, head->row == shift && head->col == shift, is_sink, len_obs, observed, &variant);
                        if (count > 0)
                        {
                            if (head->idx == GVA_NULL)
                            {
                                va_array_append(va_std_allocator, graph.nodes, ((Node2) {head->row, head->col, head->length, GVA_NULL, GVA_NULL}));
                                head->idx = va_array_length(graph.nodes) - 1;
                                fprintf(stderr, "MAKE NODE %u: (%u, %u, %u)\n", head->idx, head->row, head->col, head->length);
                                if (head->row == shift && head->col == shift)
                                {
                                    // head is the source
                                    source = *head;
                                    graph.source = source.idx;
                                    found_source = true;
                                } // if
                            } // if
                            head->incoming = MIN(head->incoming, variant.start);
                            if (variant.end + count > tail->incoming)
                            {
                                fprintf(stderr, "SPLIT %u %u %u  @ %u\n", tail->row, tail->col, tail->length, tail->incoming);
                                uint32_t const split_len = tail->incoming - tail->row;
                                va_array_append(va_std_allocator, graph.nodes,
                                                ((Node2) {tail->row, tail->col, split_len, GVA_NULL, tail->idx}));
                                split_idx = va_array_length(graph.nodes) - 1;
                                graph.nodes[tail->idx].row += split_len;
                                graph.nodes[tail->idx].col += split_len;
                                graph.nodes[tail->idx].length -= split_len;
                                tail->incoming = variant.end + count;
                                fprintf(stderr, VAR_FMT "\n", print_variant(variant, observed));
                                if (count > 1)
                                {
                                    // incoming edge
                                    fprintf(stderr, "INCOMING -> BOTH\n");
                                    va_array_append(va_std_allocator, graph.edges, ((Edge2) {split_idx, graph.nodes[head->idx].edges}));
                                    graph.nodes[head->idx].edges = va_array_length(graph.edges) - 1;
                                } // if
                                uint32_t head_head = GVA_NULL;
                                uint32_t head_tail = GVA_NULL;
                                uint32_t tail_head = GVA_NULL;
                                uint32_t tail_tail = GVA_NULL;
                                for (uint32_t i = graph.nodes[tail->idx].edges; i != GVA_NULL; i = graph.edges[i].next)
                                {
                                    // outgoing edges
                                    VA_Variant var;
                                    uint32_t const count2 = edges2(*tail, (VA_LCS_Node) {graph.nodes[graph.edges[i].tail].row, graph.nodes[graph.edges[i].tail].col, graph.nodes[graph.edges[i].tail].length, 0, 0}, tail->row == shift && tail->col == shift, graph.nodes[graph.edges[i].tail].edges == GVA_NULL && graph.nodes[graph.edges[i].tail].lambda == GVA_NULL, len_obs, observed, &var);
                                    fprintf(stderr, "    -> %u: " VAR_FMT " x %u ::\n", graph.edges[i].tail, print_variant(var, observed), count2);
                                    if (var.start > variant.end + count - 1)
                                    {
                                        fprintf(stderr, "    -> TAIL\n");
                                        if (tail_tail != GVA_NULL)
                                        {
                                            graph.edges[tail_tail].next = i;
                                        } // if
                                        if (tail_head == GVA_NULL)
                                        {
                                            tail_head = i;
                                        } // if
                                        tail_tail = i;
                                    } // if
                                    else if (var.end + count2 - 1 < variant.end + count)
                                    {
                                        fprintf(stderr, "    -> HEAD\n");
                                        if (head_tail != GVA_NULL)
                                        {
                                            graph.edges[head_tail].next = i;
                                        } // if
                                        if (head_head == GVA_NULL)
                                        {
                                            head_head = i;
                                        } // if
                                        head_tail = i;
                                    } // if
                                    else
                                    {
                                        fprintf(stderr, "    -> BOTH\n");
                                        if (tail_tail != GVA_NULL)
                                        {
                                            graph.edges[tail_tail].next = i;
                                        } // if
                                        if (tail_head == GVA_NULL)
                                        {
                                            tail_head = i;
                                        } // if
                                        tail_tail = i;
                                        va_array_append(va_std_allocator, graph.edges, ((Edge2) {graph.edges[i].tail, head_head}));
                                        head_head = va_array_length(graph.edges) - 1;
                                    } // else
                                } // for
                                if (head_tail != GVA_NULL)
                                {
                                    graph.edges[head_tail].next = GVA_NULL;
                                } // if
                                if (tail_tail != GVA_NULL)
                                {
                                    graph.edges[tail_tail].next = GVA_NULL;
                                } // if
                                graph.nodes[split_idx].edges = head_head;
                                graph.nodes[tail->idx].edges = tail_head;
                            } // if
                            va_array_append(va_std_allocator, graph.edges, ((Edge2) {tail->idx, graph.nodes[head->idx].edges}));
                            graph.nodes[head->idx].edges = va_array_length(graph.edges) - 1;
                            fprintf(stderr, "MAKE EDGE\n");
                        } // if
                        else
                        {
                           // FIXME: skip the remainder of this (head) level
                        } // else
                    } // else
                } // for h_j
                is_sink = false;
                tail->idx = split_idx;
            } // for h_i
            if (!found_source && tail->length >= len_lcs - t_i)
            {
                if (source.idx == GVA_NULL)
                {
                    va_array_append(va_std_allocator, graph.nodes, ((Node2) {source.row, source.col, source.length, GVA_NULL, GVA_NULL}));
                    source.idx = va_array_length(graph.nodes) - 1;
                    fprintf(stderr, "MAKE NODE %u: (%u, %u, %u)\n", source.idx, source.row, source.col, source.length);
                    graph.source = source.idx;
                } // if
                fprintf(stderr, "SOURCE    %u %u %u: ", source.row, source.col, source.length);
                // is this edges2 check always true?
                if (edges2(source, *tail, true, is_sink, len_obs, observed, &(VA_Variant) {0, 0, 0, 0}))
                {
                    va_array_append(va_std_allocator, graph.edges, ((Edge2) {tail->idx, graph.nodes[source.idx].edges}));
                    graph.nodes[source.idx].edges = va_array_length(graph.edges) - 1;
                    fprintf(stderr, "MAKE EDGE\n");
                } // if
            } // if
        } // for t_j
        // dealloc lcs pos level
        lcs_nodes[len_lcs - t_i - 1] = va_array_destroy(va_std_allocator, lcs_nodes[len_lcs - t_i - 1]);
    } // for t_i

    // deallocation
    lcs_nodes = va_std_allocator.alloc(va_std_allocator.context, lcs_nodes, len_lcs, 0);

    fprintf(stderr, "#nodes: %zu\n#edges: %zu\n", va_array_length(graph.nodes), va_array_length(graph.edges));
    fprintf(stderr, "source: %u\n", graph.source);
    for (size_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        fprintf(stderr, "%zu: (%u, %u, %u):\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length);
        if (graph.nodes[i].lambda != GVA_NULL)
        {
            fprintf(stderr, "    (%u, %u, %u): lambda\n", graph.nodes[graph.nodes[i].lambda].row, graph.nodes[graph.nodes[i].lambda].col, graph.nodes[graph.nodes[i].lambda].length);
        } // if
        for (uint32_t j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next)
        {
            fprintf(stderr, "    (%u, %u, %u): ", graph.nodes[graph.edges[j].tail].row, graph.nodes[graph.edges[j].tail].col, graph.nodes[graph.edges[j].tail].length);
            edges2(((VA_LCS_Node) {graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, 0, 0}),
                   ((VA_LCS_Node) {graph.nodes[graph.edges[j].tail].row, graph.nodes[graph.edges[j].tail].col, graph.nodes[graph.edges[j].tail].length, 0, 0}),
                   graph.nodes[i].row == shift && graph.nodes[i].col == shift, graph.nodes[graph.edges[j].tail].edges == GVA_NULL, len_obs, observed, &(VA_Variant) {0, 0, 0, 0});
        } // for
    } // for

    //bfs_traversal2(graph, len_obs, observed);

    return graph;
} // build


void
to_json2(Graph2 const graph, size_t const len_obs, char const observed[static len_obs])
{
    printf("{\n    \"source\": \"s%u\",\n    \"nodes\": {\n", graph.source);
    for (uint32_t i = 0; i < va_array_length(graph.nodes); ++i)
    {
        printf("        \"s%u\": {\"row\": %u, \"col\": %u, \"length\": %u}%s\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, i < va_array_length(graph.nodes) - 1 ? "," : "");
    } // for
    printf("    },\n");
    bfs_traversal2(graph, len_obs, observed);
    printf("}");
} // to_json2
