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


static inline void
split(VA_LCS_Node2 const* const node,
      uint32_t const in_count, VA_Variant const incoming,
      uint32_t const split_point,
      Graph2* const graph)
{
/*
    fprintf(stderr, "\n***SPLIT (%u, %u, %u) because of first outgoing edge at %u\n", node->row, node->col, node->length, node->outgoing);
    fprintf(stderr, "fsp: %u lsp: %u\n", MIN(split_point - 1, incoming.end + in_count - 1), MAX(incoming.end, node->outgoing));
*/
    uint32_t pos = 0;
    for (pos = MIN(split_point - 1, incoming.end + in_count - 1); pos >= MAX(incoming.end, node->outgoing); --pos)
    {
        uint32_t head_head = GVA_NULL;
        uint32_t head_tail = GVA_NULL;
        uint32_t tail_head = GVA_NULL;
        uint32_t tail_tail = GVA_NULL;
        // FIXME: only once after creating all the nodes
        for (uint32_t j = graph->nodes[node->idx].edges; j != GVA_NULL; j = graph->edges[j].next)
        {
            Node2 const* const tail = &graph->nodes[graph->edges[j].tail];
            VA_Variant outgoing;
            uint32_t const out_count = edges2(
                (VA_LCS_Node) {node->row, node->col, node->length, -1, GVA_NULL, 0},
                (VA_LCS_Node) {tail->row, tail->col, tail->length, -1, GVA_NULL, 0}, false, tail->edges == GVA_NULL && tail->lambda == GVA_NULL, &outgoing);
            //fprintf(stderr, "    %u: (%u, %u, %u): x %u  ", graph->edges[j].tail, tail->row, tail->col, tail->length, out_count);

            // TODO: can we avoid reordering the tail list? Yes? XOR-linked?
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

                //fprintf(stderr, "tail\n");
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

                //fprintf(stderr, "head\n");
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
                if (head_tail == GVA_NULL)
                {
                    head_tail = head_head;
                } // if

                //fprintf(stderr, "both\n");
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
        graph->nodes[node->idx].edges = head_head;
        uint32_t const length = pos - node->row;

        va_array_append(va_std_allocator, graph->nodes, ((Node2) {node->row + length, node->col + length, graph->nodes[node->idx].length - length, tail_head, graph->nodes[node->idx].lambda}));
        graph->nodes[node->idx].lambda = va_array_length(graph->nodes) - 1;
        graph->nodes[node->idx].length = length;

    /*
        Node2 const* const lambda = &graph->nodes[graph->nodes[node->idx].lambda];
        fprintf(stderr, "(%u, %u, %u) l> (%u, %u, %u)  @ %u\n",
                        graph->nodes[node->idx].row, graph->nodes[node->idx].col, graph->nodes[node->idx].length,
                        lambda->row, lambda->col, lambda->length, pos);
                 */
    } // for

    //fprintf(stderr, "***\n\n");
} // split


static inline void
distribute(VA_LCS_Node2 const* const tail, VA_LCS_Node2 const* const head,
           uint32_t const count, VA_Variant const variant,
           Graph2* const graph)
{
    uint32_t lambda = tail->idx;
    for (uint32_t i = 0; i < count; ++i)
    {
        //fprintf(stderr, "LAMBDA: %u\n", lambda);
        //fprintf(stderr, "    variant: %u\n", variant.end + i);
        bool is_lambda = false;
        while (graph->nodes[lambda].lambda != GVA_NULL &&
               variant.end + i >= graph->nodes[lambda].row + graph->nodes[lambda].length)
        {
            lambda = graph->nodes[lambda].lambda;
            is_lambda = true;
        } // while

        if (is_lambda || i == 0)
        {
            va_array_append(va_std_allocator, graph->edges, ((Edge2) {lambda, graph->nodes[head->idx].edges}));
            graph->nodes[head->idx].edges = va_array_length(graph->edges) - 1;
        } // if

        //fprintf(stderr, "    next: %d\n", graph->nodes[lambda].lambda);
        if (graph->nodes[lambda].lambda == GVA_NULL || variant.end + count <= graph->nodes[lambda].row + graph->nodes[lambda].length)
        {
            break;
        } // if
    } // for
} // distribute


static inline uint32_t
edges3(VA_LCS_Node2 const* const head, VA_LCS_Node2 const* const tail,
       bool const is_source, bool const is_sink,
       VA_Variant* const variant)
{
    uint32_t const h_min = head->lcs_pos - (head->length - 1);
    uint32_t const t_min = tail->lcs_pos - (tail->length - 1);

    uint32_t const start = MAX(h_min + 1 - is_source, t_min);
    uint32_t const end = MIN(head->lcs_pos + 1, tail->lcs_pos + is_sink);

    if (start > end)
    {
        return 0;
    } // if

    uint32_t const h_offset = (head->length - 1) - (head->lcs_pos - start);
    uint32_t const t_offset = (tail->length - 1) - (tail->lcs_pos - start);

    if (head->row + h_offset > tail->row + t_offset || head->col + h_offset > tail->col + t_offset)
    {
        return 0;
    } // if

    *variant = (VA_Variant) {head->row + h_offset, tail->row + t_offset, head->col + h_offset, tail->col + t_offset};
    return end - start + 1;
} // edges3


Graph2
build(size_t const len_ref, char const reference[static restrict len_ref],
      size_t const len_obs, char const observed[static restrict len_obs],
      size_t const shift)
{
    uint32_t last_idx = -1;
    VA_LCS_Node2* lcs_nodes = va_edit2(va_std_allocator, len_ref, reference, len_obs, observed, &last_idx);

    Graph2 graph = {NULL, NULL, GVA_NULL};

    if (lcs_nodes == NULL)
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

/*
    for (size_t i = 0; i < va_array_length(lcs_nodes); ++i)
    {
        fprintf(stderr, "%2zu: (%u, %u, %u)  %u  %2d  %2d\n", i, lcs_nodes[i].row, lcs_nodes[i].col, lcs_nodes[i].length, lcs_nodes[i].lcs_pos, lcs_nodes[i].next, lcs_nodes[i].inext);
    } // for
    fprintf(stderr, "%u\n", last_idx);
*/
    VA_LCS_Node2 source = {shift, shift, 0, -1, GVA_NULL, GVA_NULL, GVA_NULL, -1};
    bool found_source = false;
    size_t const len = va_array_length(lcs_nodes);
    VA_LCS_Node2* const last = &lcs_nodes[last_idx];
    if (last->row + last->length == len_ref + shift && last->col + last->length == len_obs + shift)
    {
        va_array_append(va_std_allocator, graph.nodes, ((Node2) {last->row, last->col, last->length, GVA_NULL, GVA_NULL}));
        last->idx = va_array_length(graph.nodes) - 1;
        if (last->row == shift && last->col == shift)
        {
            found_source = true;
            source.idx = last->idx;
        } // if
    } // if
    else
    {
        VA_LCS_Node2 const sink = {len_ref, len_obs, 0, lcs_nodes[last_idx].lcs_pos, GVA_NULL, GVA_NULL, GVA_NULL, -1};
        va_array_append(va_std_allocator, graph.nodes, ((Node2) {sink.row, sink.col, sink.length, GVA_NULL, GVA_NULL}));
        uint32_t const sink_idx = va_array_length(graph.nodes) - 1;
        for (uint32_t j = last_idx; j != GVA_NULL && lcs_nodes[j].lcs_pos == lcs_nodes[last_idx].lcs_pos; j = lcs_nodes[j].next)
        {
            VA_LCS_Node2* const head = &lcs_nodes[j];
            va_array_append(va_std_allocator, graph.edges, ((Edge2) {sink_idx, GVA_NULL}));
            va_array_append(va_std_allocator, graph.nodes, ((Node2) {head->row, head->col, head->length, va_array_length(graph.edges) - 1, GVA_NULL}));
            head->idx = va_array_length(graph.nodes) - 1;
            head->outgoing = head->row + head->length;
            if (head->row == shift && head->col == shift)
            {
                found_source = true;
                source.idx = head->idx;
                source.outgoing = head->outgoing;
            } // if

/*
            VA_Variant variant;
            uint32_t const count = edges3(head, &sink, head->row == shift && head->col == shift, true, &variant);
            fprintf(stderr, "(%u, %u, %u) -> (%u, %u, %u)  " VAR_FMT " x %u\n", head->row, head->col, head->length, sink.row, sink.col, sink.length, print_variant(variant, observed), count);
*/
        } // for
    } // else
//    fprintf(stderr, "\n");

    for (size_t i = 0; i < len; ++i)
    {
        VA_LCS_Node2 const* const restrict tail = &lcs_nodes[len - i - 1];
        uint32_t const idx = len - i - 1;

        if (lcs_nodes[idx].inext != GVA_NULL)
        {
            lcs_nodes[lcs_nodes[idx].inext].next = lcs_nodes[idx].next;
        } // if
        else
        {
            last_idx = lcs_nodes[idx].next;
        } // else
        if (lcs_nodes[idx].next != GVA_NULL)
        {
            lcs_nodes[lcs_nodes[idx].next].inext = lcs_nodes[idx].inext;
        } // if

        if (tail->idx == GVA_NULL)
        {
            continue;
        } // if
        bool const is_sink = tail->row + tail->length == len_ref + shift && tail->col + tail->length == len_obs + shift;

        uint32_t split_point = -1;

        uint32_t const lcs_low = tail->length > tail->lcs_pos ? 0 : tail->lcs_pos - tail->length;
        //fprintf(stderr, "(%u, %u, %u)  @ %u  low: %u\n", tail->row, tail->col, tail->length, tail->lcs_pos, lcs_low);

        for (uint32_t j = last_idx; j != GVA_NULL && lcs_nodes[j].lcs_pos >= lcs_low; j = lcs_nodes[j].next)
        {
            //fprintf(stderr, "  %u: (%u, %u, %u)  @ %u\n", j, lcs_nodes[j].row, lcs_nodes[j].col, lcs_nodes[j].length, lcs_nodes[j].lcs_pos);
            VA_LCS_Node2* const head = &lcs_nodes[j];
            bool const is_source = head->row == shift && head->col == shift;

            VA_Variant variant;
            uint32_t const count = edges3(head, tail, is_source, is_sink, &variant);
            if (count == 0)
            {
                continue;
            } // if

            if (head->idx == GVA_NULL)
            {
                va_array_append(va_std_allocator, graph.nodes, ((Node2) {head->row, head->col, head->length, GVA_NULL, GVA_NULL}));
                head->idx = va_array_length(graph.nodes) - 1;
                if (is_source)
                {
                    found_source = true;
                    source.idx = head->idx;
                } // if
            } // if

            //fprintf(stderr, "(%u, %u, %u) -> (%u, %u, %u)  " VAR_FMT " x %u\n", head->row, head->col, head->length, tail->row, tail->col, tail->length, print_variant(variant, observed), count);

            if (variant.end + count > tail->outgoing && split_point != tail->outgoing)
            {
                split(tail, count, variant, split_point, &graph);
                split_point = MIN(split_point, MAX(variant.end, tail->outgoing));
            } // if
            distribute(tail, head, count, variant, &graph);
            head->outgoing = MIN(head->outgoing, variant.start);
            source.outgoing = MIN(source.outgoing, variant.start);
        } // for

        if (!found_source && tail->length > tail->lcs_pos)
        {
            if (source.idx == GVA_NULL)
            {
                va_array_append(va_std_allocator, graph.nodes, ((Node2) {source.row, source.col, source.length, GVA_NULL, GVA_NULL}));
                source.idx = va_array_length(graph.nodes) - 1;
                source.outgoing = 0;
            } // if
            va_array_append(va_std_allocator, graph.edges, ((Edge2) {tail->idx, graph.nodes[source.idx].edges}));
            graph.nodes[source.idx].edges = va_array_length(graph.edges) - 1;

/*
            VA_Variant variant;
            uint32_t const count = edges3(&source, tail, true, is_sink, &variant);
            fprintf(stderr, "(%u, %u, %u) -> (%u, %u, %u)  " VAR_FMT " x %u\n", source.row, source.col, source.length, tail->row, tail->col, tail->length, print_variant(variant, observed), count);
*/
        } // if
        //fprintf(stderr, "\n");

    } // for

    lcs_nodes = va_array_destroy(va_std_allocator, lcs_nodes);

    //fprintf(stderr, "SOURCE first outgoing edge at %d\n", source.outgoing);
    if (source.outgoing != (uint32_t) -1)
    {
        graph.nodes[source.idx].row += source.outgoing;
        graph.nodes[source.idx].col += source.outgoing;
        graph.nodes[source.idx].length -= source.outgoing;
    } // if
    else
    {
        graph.nodes[source.idx].length = 0;
    } // else

    graph.source = source.idx;
    //print_graph(graph, len_obs, observed);

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
