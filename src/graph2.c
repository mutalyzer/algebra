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
                uint32_t const count2 = edges2(((VA_LCS_Node) {graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, 0, 0}),
                                               ((VA_LCS_Node) {graph.nodes[graph.edges[j].tail].row, graph.nodes[graph.edges[j].tail].col, graph.nodes[graph.edges[j].tail].length, 0, 0}),
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
                    ((VA_LCS_Node) {graph.nodes[i].row, graph.nodes[i].col, graph.nodes[i].length, 0, 0}),
                    ((VA_LCS_Node) {graph.nodes[graph.edges[j].tail].row, graph.nodes[graph.edges[j].tail].col,
                                    graph.nodes[graph.edges[j].tail].length, 0, 0}),
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
            uint32_t const out_count = edges2(*node, (VA_LCS_Node) {tail->row, tail->col, tail->length, -1, GVA_NULL}, false, tail->edges == GVA_NULL && tail->lambda == GVA_NULL, &outgoing);
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
    VA_LCS_Node** lcs_nodes = NULL;
    size_t const len_lcs = va_edit(va_std_allocator, len_ref, reference, len_obs, observed, &lcs_nodes);

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

    VA_LCS_Node source = {shift, shift, 0, -1, GVA_NULL};
    bool found_source = false;
    VA_LCS_Node* const last = &lcs_nodes[len_lcs - 1][va_array_length(lcs_nodes[len_lcs - 1]) - 1];
    if (last->row + last->length != len_ref + shift || last->col + last->length != len_obs + shift)
    {
        VA_LCS_Node const sink = {len_ref, len_obs, 0, -1, GVA_NULL};
        va_array_append(va_std_allocator, graph.nodes, ((Node2) {sink.row, sink.col, sink.length, GVA_NULL, GVA_NULL}));
        uint32_t const sink_idx = va_array_length(graph.nodes) - 1;
        size_t const h_len = va_array_length(lcs_nodes[len_lcs - 1]);
        for (size_t i = 0; i < h_len; ++i)
        {
            VA_LCS_Node* const head = &lcs_nodes[len_lcs - 1][h_len - i - 1];
            if (head->row == shift && head->col == shift)
            {
                found_source = true;
            } // if

            va_array_append(va_std_allocator, graph.edges, ((Edge2) {sink_idx, GVA_NULL}));
            va_array_append(va_std_allocator, graph.nodes, ((Node2) {head->row, head->col, head->length, va_array_length(graph.edges) - 1, GVA_NULL}));
            head->idx = va_array_length(graph.nodes) - 1;
            head->incoming = head->row + head->length;

            VA_Variant variant;
            uint32_t const count = edges2(*head, sink, head->row == shift && head->col == shift, true, &variant);
            fprintf(stderr, "SINK (%u, %u, %u) -> (%u, %u, %u)  " VAR_FMT " x %u\n", head->row, head->col, head->length, sink.row, sink.col, sink.length, print_variant(variant, observed), count);
        } // for
    } // if
    else
    {
        va_array_append(va_std_allocator, graph.nodes, ((Node2) {last->row, last->col, last->length, GVA_NULL, GVA_NULL}));
        last->idx = va_array_length(graph.nodes) - 1;

        if (last->row == shift && last->col == shift)
        {
            found_source = true;
            source.idx = last->idx;
            last->incoming = 0;
            source.incoming = 0;
        } // if
    } // else

    size_t* x = va_array_init(va_std_allocator, len_lcs, sizeof(*x));
    for (size_t i = 0; i < len_lcs; ++i)
    {
        x[i] = va_array_length(lcs_nodes[i]);
    } // for

    for (size_t t_i = 0; t_i < len_lcs; ++t_i)
    {
        size_t const t_len = va_array_length(lcs_nodes[len_lcs - t_i - 1]);
        for (size_t t_j = 0; t_j < t_len; ++t_j)
        {
            VA_LCS_Node* const tail = &lcs_nodes[len_lcs - t_i - 1][t_len - t_j - 1];
            if (tail->idx == GVA_NULL)
            {
                continue;
            } // if

            if (tail->row == shift && tail->col == shift)
            {
                source.idx = tail->idx;
                source.incoming = tail->incoming;
            } // if

            bool const is_sink = tail->row + tail->length == len_ref + shift && tail->col + tail->length == len_obs + shift;
            size_t const len = MIN(len_lcs - t_i, tail->length + 1);
            bool converse_stop = false;
            size_t append_level = 0;
            for (size_t h_i = 0; !converse_stop && h_i < len; ++h_i)
            {
                size_t const h_len = MIN(h_i == 0 ? t_len - t_j - 1: va_array_length(lcs_nodes[len_lcs - t_i - h_i - 1]), x[len_lcs - t_i - h_i - 1]);
                for (size_t h_j = 0; h_j < h_len; ++h_j)
                {
                    VA_LCS_Node* const head = &lcs_nodes[len_lcs - t_i - h_i - 1][h_len - h_j - 1];
                    if (head->row == shift && head->col == shift)
                    {
                        found_source = true;
                    } // if

                    fprintf(stderr, "(%u, %u, %u) vs (%u, %u, %u)\n", head->row, head->col, head->length, tail->row, tail->col, tail->length);

                    if (head->row + head->length + h_i <= tail->row + tail->length &&
                        head->col + head->length + h_i <= tail->col + tail->length)
                    {
                        VA_Variant variant;
                        uint32_t const count = edges2(*head, *tail, head->row == shift && head->col == shift, is_sink, &variant);
                        if (count > 0)
                        {
                            if (head->idx == GVA_NULL)
                            {
                                va_array_append(va_std_allocator, graph.nodes, ((Node2) {head->row, head->col, head->length, GVA_NULL, GVA_NULL}));
                                head->idx = va_array_length(graph.nodes) - 1;
                            } // if

                            if (variant.end + count > tail->incoming)
                            {
                                split(tail, head, count, variant, &graph);
                            } // if
                            else
                            {
                                uint32_t lambda = tail->idx;
                                for (uint32_t i = 0; i < count; ++i)
                                {
                                    fprintf(stderr, "LAMBDA: %u\n", lambda);
                                    fprintf(stderr, "    variant: %u\n", variant.end + i);
                                    bool is_lambda = false;
                                    while (graph.nodes[lambda].lambda != GVA_NULL &&
                                           variant.end + i >= graph.nodes[lambda].row + graph.nodes[lambda].length)
                                    {
                                        lambda = graph.nodes[lambda].lambda;
                                        is_lambda = true;
                                    } // if

                                    if (i == 0 || is_lambda)
                                    {
                                        va_array_append(va_std_allocator, graph.edges, ((Edge2) {lambda, graph.nodes[head->idx].edges}));
                                        graph.nodes[head->idx].edges = va_array_length(graph.edges) - 1;
                                    } // if

                                    fprintf(stderr, "    next: %u\n", graph.nodes[lambda].lambda);
                                    if (graph.nodes[lambda].lambda == GVA_NULL || variant.end + count <= graph.nodes[lambda].row + graph.nodes[lambda].length)
                                    {
                                        break;
                                    } // if
                                } // for
                                head->incoming = MIN(head->incoming, variant.start);
                            } // else

                            fprintf(stderr, "(%u, %u, %u) -> (%u, %u, %u)  " VAR_FMT " x %u\n", head->row, head->col, head->length, tail->row, tail->col, tail->length, print_variant(variant, observed), count);
                        } // if
                    } // if
                    else if (tail->row + tail->length <= head->row + head->length + h_i &&
                             tail->col + tail->length <= head->col + head->length + h_i)
                    {
                        VA_Variant variant;
                        uint32_t const count = edges2(*tail, *head, tail->row == shift && tail->col == shift, false, &variant);
                        if (count > 0)
                        {
                            if (head->idx == GVA_NULL)
                            {
                                fprintf(stderr, "(%u, %u, %u) not in the graph yet\n", head->row, head->col, head->length);
                                // FIXME: why is this special?
                                if (h_i > 0)
                                {
                                    fprintf(stderr, "Marking to stop\n");
                                    append_level = h_i;
                                    converse_stop = true;
                                } // if
                            } // if
                            else if (variant.end + count > head->incoming)
                            {
                                split(head, tail, count, variant, &graph);
                            } // if
                            else
                            {
                                // FIXME: add to lambda chain?
                                va_array_append(va_std_allocator, graph.edges, ((Edge2) {head->idx, graph.nodes[tail->idx].edges}));
                                graph.nodes[tail->idx].edges = va_array_length(graph.edges) - 1;
                                tail->incoming = MIN(tail->incoming, variant.start);
                                if (source.idx == tail->idx)
                                {
                                    source.incoming = tail->incoming;
                                } // if
                            } // else

                            fprintf(stderr, "CONVERSE (%u, %u, %u) -> (%u, %u, %u)  " VAR_FMT " x %u\n", tail->row, tail->col, tail->length, head->row, head->col, head->length, print_variant(variant, observed), count);
                        } // if
                    } // if

                } // for h_j
            } // for h_i
            if (converse_stop)
            {
                fprintf(stderr, "Append (%u, %u, %zu) at level: %zu\n", tail->row, tail->col, tail->length - append_level, len_lcs - t_i - append_level - 1);
                va_array_append(va_std_allocator, lcs_nodes[len_lcs - t_i - append_level - 1], ((VA_LCS_Node) {tail->row, tail->col, tail->length - append_level, tail->incoming, tail->idx}));
            } // if
            else if (!found_source && tail->length >= len_lcs - t_i)
            {
                if (source.idx == GVA_NULL)
                {
                    va_array_append(va_std_allocator, graph.nodes, ((Node2) {source.row, source.col, source.length, GVA_NULL, GVA_NULL}));
                    source.idx = va_array_length(graph.nodes) - 1;
                } // if
                va_array_append(va_std_allocator, graph.edges, ((Edge2) {tail->idx, graph.nodes[source.idx].edges}));
                graph.nodes[source.idx].edges = va_array_length(graph.edges) - 1;
                source.incoming = MIN(source.incoming, source.row + source.length);

                VA_Variant variant;
                uint32_t const count = edges2(source, *tail, true, is_sink, &variant);
                fprintf(stderr, "SOURCE (%u, %u, %u) -> (%u, %u, %u)  " VAR_FMT " x %u\n", source.row, source.col, source.length, tail->row, tail->col, tail->length, print_variant(variant, observed), count);
            } // if
        } // for t_j
        lcs_nodes[len_lcs - t_i - 1] = va_array_destroy(va_std_allocator, lcs_nodes[len_lcs - t_i - 1]);
    } // for t_i
    lcs_nodes = va_std_allocator.alloc(va_std_allocator.context, lcs_nodes, len_lcs, 0);

    x = va_array_destroy(va_std_allocator, x);

    fprintf(stderr, "SOURCE first outgoing edge at %u\n", source.incoming);
    graph.nodes[source.idx].row += source.incoming;
    graph.nodes[source.idx].col += source.incoming;
    graph.nodes[source.idx].length -= source.incoming;

    graph.source = source.idx;

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
