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
       VA_Variant* const variant)
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
        //fprintf(stderr, "NO EDGE\n");
        return 0;
    } // if

    uint32_t const count = MIN(head_length - head_offset, tail_length - tail_offset - 1) + 1;

    *variant = (VA_Variant) {row + head_offset, tail.row + tail_offset, col + head_offset, tail.col + tail_offset};

    //fprintf(stderr, VAR_FMT " x %u\n", print_variant((*variant), observed), count);
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
                uint32_t const count2 = edges2(((VA_LCS_Node) {graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, 0, 0, 0}),
                                               ((VA_LCS_Node) {graph.nodes[graph.edges[j].tail].row, graph.nodes[graph.edges[j].tail].col, graph.nodes[graph.edges[j].tail].length, 0, 0, 0}),
                                               graph.nodes[head].row == graph.nodes[graph.source].row && graph.nodes[head].col == graph.nodes[graph.source].col, graph.nodes[graph.edges[j].tail].edges == GVA_NULL, &variant);

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


static void
print_graph(Graph2 const graph, size_t const len_obs, char const observed[static len_obs]) {
    fprintf(stderr, "#nodes: %zu\n#edges: %zu\n", va_array_length(graph.nodes), va_array_length(graph.edges));
    fprintf(stderr, "source: %u\n", graph.source);
    if (graph.nodes == NULL)
    {
        return;
    } // if
    for (size_t i = 0; i < va_array_length(graph.nodes); ++i) {
        fprintf(stderr, "%zu: (%u, %u, %u):\n", i, graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length);
        if (graph.nodes[i].lambda != GVA_NULL) {
            fprintf(stderr, "    (%u, %u, %u): lambda\n", graph.nodes[graph.nodes[i].lambda].row, graph.nodes[graph.nodes[i].lambda].col, graph.nodes[graph.nodes[i].lambda].length);
        } // if
        for (uint32_t j = graph.nodes[i].edges; j != GVA_NULL; j = graph.edges[j].next) {
            fprintf(stderr, "    %u: (%u, %u, %u): ", graph.edges[j].tail, graph.nodes[graph.edges[j].tail].row, graph.nodes[graph.edges[j].tail].col, graph.nodes[graph.edges[j].tail].length);
            VA_Variant variant;
            uint32_t const count = edges2(
                    ((VA_LCS_Node) {graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, 0, 0, 0}),
                    ((VA_LCS_Node) {graph.nodes[graph.edges[j].tail].row, graph.nodes[graph.edges[j].tail].col,
                                    graph.nodes[graph.edges[j].tail].length, 0, 0, 0}),
                    i == graph.source,
                    graph.nodes[graph.edges[j].tail].edges == GVA_NULL && graph.nodes[graph.edges[j].tail].edges == GVA_NULL,
                    &variant);
            fprintf(stderr, VAR_FMT " x %u\n", print_variant(variant, observed), count);
        } // for
    } // for
} // print_graph


static void
split(VA_LCS_Node* const node, VA_LCS_Node* const head, uint32_t const in_count, VA_Variant const incoming, Graph2* const graph)
{
    fprintf(stderr, "\n***SPLIT (%u, %u, %u) because of first outgoing edge at %u\n", node->row, node->col, node->length, node->incoming);

    uint32_t pos = 0;
    for (uint32_t i = 0; i < incoming.end + in_count - node->incoming; ++i)
    {
        pos = incoming.end + in_count - i - 1;

        uint32_t head_head = GVA_NULL;
        uint32_t head_tail = GVA_NULL;
        uint32_t tail_head = GVA_NULL;
        uint32_t tail_tail = GVA_NULL;
        for (uint32_t j = graph->nodes[node->idx].edges; j != GVA_NULL; j = graph->edges[j].next)
        {
            Node2 const* const tail = &graph->nodes[graph->edges[j].tail];
            VA_Variant outgoing;
            uint32_t const out_count = edges2(*node, (VA_LCS_Node) {tail->row, tail->col, tail->length, -1, GVA_NULL, 0}, false, tail->edges == GVA_NULL && tail->lambda == GVA_NULL, &outgoing);
            fprintf(stderr, "    %u: (%u, %u, %u): x %u  ", graph->edges[j].tail, tail->row, tail->col, tail->length, out_count);

            // TODO: can we avoid reordering the tail list?
            if (outgoing.start > pos)
            {
                if (tail_tail != GVA_NULL)
                {
                    graph->edges[tail_tail].next = j;
                } // if
                if (tail_head == GVA_NULL)
                {
                    tail_head = j;
                } // if
                tail_tail = j;

                fprintf(stderr, "tail\n");
            } // if
            else if (outgoing.start + out_count - 1 <= pos)
            {
                if (head_tail != GVA_NULL)
                {
                    graph->edges[head_tail].next = j;
                } // if
                if (head_head == GVA_NULL)
                {
                    head_head = j;
                } // if
                head_tail = j;

                fprintf(stderr, "head\n");
            } // if
            else
            {
                if (tail_tail != GVA_NULL)
                {
                    graph->edges[tail_tail].next = j;
                } // if
                if (tail_head == GVA_NULL)
                {
                    tail_head = j;
                } // if
                tail_tail = j;
                va_array_append(va_std_allocator, graph->edges, ((Edge2) {graph->edges[j].tail, head_head}));
                head_head = va_array_length(graph->edges) - 1;

                fprintf(stderr, "both\n");
            } // else
        } // for
        if (head_tail != GVA_NULL)
        {
            graph->edges[head_tail].next = GVA_NULL;
        } // if
        if (tail_tail != GVA_NULL)
        {
            graph->edges[tail_tail].next = GVA_NULL;
        } // if
        graph->nodes[node->idx].edges = tail_head;

        graph->nodes[node->idx].row += (pos - node->row);
        graph->nodes[node->idx].col += (pos - node->row);
        graph->nodes[node->idx].length -= (pos - node->row);

        va_array_append(va_std_allocator, graph->edges, ((Edge2) {node->idx, graph->nodes[head->idx].edges}));
        graph->nodes[head->idx].edges = va_array_length(graph->edges) - 1;

        va_array_append(va_std_allocator, graph->nodes, ((Node2) {node->row, node->col, pos - node->row, head_head, node->idx}));
        node->idx = va_array_length(graph->nodes) - 1;

        Node2 const* const lambda = &graph->nodes[graph->nodes[node->idx].lambda];
        fprintf(stderr, "(%u, %u, %u) l> (%u, %u, %u)  @ %u\n", graph->nodes[node->idx].row, graph->nodes[node->idx].col, graph->nodes[node->idx].length,
                                                          lambda->row, lambda->col, lambda->length, pos);
    } // for

    if (incoming.end < pos)
    {
        va_array_append(va_std_allocator, graph->edges, ((Edge2) {node->idx, graph->nodes[head->idx].edges}));
        graph->nodes[head->idx].edges = va_array_length(graph->edges) - 1;
    } // if
    node->incoming = -1;  // FIXME: real value?
    head->incoming = incoming.start;

    fprintf(stderr, "***\n\n");
} // split


// FIXME: incoming means the first *outgoing* edge
Graph2
build(size_t const len_ref, char const reference[static len_ref],
      size_t const len_obs, char const observed[static len_obs],
      size_t const shift)
{
    VA_LCS_Node2* lcs_nodes = NULL;
    uint32_t* lcs_index = NULL;
    size_t const len_lcs = va_edit2(va_std_allocator, len_ref, reference, len_obs, observed, &lcs_nodes, &lcs_index);

    Graph2 graph = {NULL, NULL, GVA_NULL};

    if (len_lcs == 0 || lcs_nodes == NULL)
    {
        va_array_append(va_std_allocator, graph.nodes, ((Node2) {len_ref, len_obs, 0, GVA_NULL, GVA_NULL}));
        uint32_t const sink = va_array_length(graph.nodes) - 1;
        if (len_ref == 0 && len_obs == 0)
        {
            graph.source = sink;
            return graph;
        } // if

        va_array_append(va_std_allocator, graph.edges, ((Edge2) {sink, GVA_NULL}));
        va_array_append(va_std_allocator, graph.nodes, ((Node2) {shift, shift, 0, va_array_length(graph.edges) - 1, GVA_NULL}));
        graph.source = va_array_length(graph.nodes) - 1;
        return graph;
    } // if

    for (size_t i = 0; i < va_array_length(lcs_nodes); ++i)
    {
        fprintf(stderr, "%zu: (%u, %u, %u), %d\n", i, lcs_nodes[i].row, lcs_nodes[i].col, lcs_nodes[i].length, lcs_nodes[i].prev);
    } // for

    for (size_t i = 0; i < len_lcs; ++i)
    {
        fprintf(stderr, "%zu: %d\n", i, lcs_index[i]);
    } // for

    VA_LCS_Node2* const last = &lcs_nodes[lcs_index[len_lcs - 1]];
    if (last->row + last->length != len_ref + shift || last->col + last->length != len_obs + shift)
    {
        fprintf(stderr, "Nieuwe sink\n");
        VA_LCS_Node2 const sink = {len_ref, len_obs, 0, len_lcs, GVA_NULL, GVA_NULL};
        va_array_append(va_std_allocator, graph.nodes, ((Node2) {sink.row, sink.col, sink.length, GVA_NULL, GVA_NULL}));
        uint32_t const sink_idx = va_array_length(graph.nodes) - 1;
        for (uint32_t idx = lcs_index[len_lcs - 1]; idx != GVA_NULL; idx = lcs_nodes[idx].prev)
        {
            VA_LCS_Node2* const head = &lcs_nodes[idx];
            va_array_append(va_std_allocator, graph.edges, ((Edge2) {sink_idx, GVA_NULL}));
            va_array_append(va_std_allocator, graph.nodes, ((Node2) {head->row, head->col, head->length, va_array_length(graph.edges) - 1, GVA_NULL}));
            head->idx = va_array_length(graph.nodes) - 1;
        } // for
    } // if
    else
    {
        fprintf(stderr, "Oude sink\n");
        va_array_append(va_std_allocator, graph.nodes, ((Node2) {last->row, last->col, last->length, GVA_NULL, GVA_NULL}));
        last->idx = va_array_length(graph.nodes) - 1;
    } // else

    size_t const t_len = va_array_length(lcs_nodes);
    for (size_t i = 0; i < t_len; ++i)
    {
        VA_LCS_Node2* const tail = &lcs_nodes[t_len - i - 1];
        fprintf(stderr, "(%u, %u, %u) @ %u\n", tail->row, tail->col, tail->length, tail->lcs_pos);
        if (tail->idx == GVA_NULL)
        {
            fprintf(stderr, "    Not in graph\n");
            tail->lcs_pos = GVA_NULL;
            continue;
        } // if

        bool const is_sink = tail->row + tail->length == len_ref + shift && tail->col + tail->length == len_obs + shift;
        uint32_t const h_len = MIN(tail->length, tail->lcs_pos);
        for (uint32_t j = 0; j <= h_len; ++j)
        {
            for (uint32_t k = j == 0 ? tail->prev : lcs_index[tail->lcs_pos - j]; k != GVA_NULL; k = lcs_nodes[k].prev)
            {
                VA_LCS_Node2* const head = &lcs_nodes[k];
                if (head->lcs_pos == GVA_NULL)
                {
                    continue;
                } // if

                fprintf(stderr, "    (%u, %u, %u)\n", head->row, head->col, head->length);
                if (head->row + head->length + j <= tail->row + tail->length &&
                    head->col + head->length + j <= tail->col + tail->length)
                {
                    VA_Variant variant;
                    uint32_t const count = edges2((VA_LCS_Node) {head->row, head->col, head->length}, (VA_LCS_Node) {tail->row, tail->col, tail->length}, head->row == shift && head->col == shift, is_sink, &variant);
                    if (count > 0)
                    {
                        fprintf(stderr, "      BINGO!\n");
                        if (head->idx == GVA_NULL)
                        {
                            va_array_append(va_std_allocator, graph.nodes, ((Node2) {head->row, head->col, head->length, GVA_NULL, GVA_NULL}));
                            head->idx = va_array_length(graph.nodes) - 1;
                        } // if
                        va_array_append(va_std_allocator, graph.edges, ((Edge2) {tail->idx, graph.nodes[head->idx].edges}));
                        graph.nodes[head->idx].edges = va_array_length(graph.edges);
                    } // if
                } // if
            } // for
        } // for

        if (tail->length > tail->lcs_pos)
        {
            fprintf(stderr, "    SOURCE\n");
        } // if
    } // for


    lcs_nodes = va_array_destroy(va_std_allocator, lcs_nodes);
    va_std_allocator.alloc(va_std_allocator.context, lcs_index, len_lcs * sizeof(*lcs_index), 0);

    print_graph(graph, len_obs, observed);

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
