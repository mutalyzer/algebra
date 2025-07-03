#define _XOPEN_SOURCE


#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/faststabber.h"     // Node, Stack


#define TIC(X) clock_t const X = clock()
#define TOC(X) fprintf(stderr, "time %s: %.3f s\n", (#X), (double)(clock() - (X)) / CLOCKS_PER_SEC)


Node*
stack_top(Stack const* const restrict stack)
{
    assert(stack != NULL);
    return stack->node;
} // stack_top


Stack*
stack_push(Stack* const restrict stack, Node* const restrict node)
{
    Stack* const restrict top = malloc(sizeof(*top));
    assert(top != NULL);
    top->node = node;
    top->next = stack;
    return top;
} // stack_push


Stack*
stack_pop(Stack* const restrict stack)
{
    assert(stack != NULL);
    Stack* const restrict top = stack->next;
    free(stack);
    return top;
} // stack_pop


Node*
create_node(uint32_t const start, uint32_t const end)
{
    Node* const restrict node = malloc(sizeof(*node));
    assert(node != NULL);

    node->start = start;
    node->end = end;
    node->parent = NULL;
    node->leftsibling =  NULL;
    node->rightchild = NULL;
    return node;
} // create_node


Node*
add_node(Node* const restrict last_node, Node* const restrict node)
{
    assert(node != NULL);

    // get parent
    node->parent = last_node;
    while (node->parent != NULL && node->parent->end < node->end)
    {
        node->parent = node->parent->parent;
    } // while

    assert(node->parent != NULL);

    // push new sibling
    node->leftsibling = node->parent->rightchild;
    node->parent->rightchild = node;

    return node;
} // add_node


size_t
draw_edge(Node const* const restrict node, void* const restrict arg)
{
    assert(node != NULL);
    if (node->parent != NULL)
    {
        return fprintf(arg, "n%p -> n%p[dir=none];\n", (void*) node->parent, (void*) node);
    } // if
    return 0;
} // draw_edge


size_t
draw_node(Node const* const restrict node, void* const restrict arg)
{
    return fprintf(arg, "n%p[label=\"%u&ndash;%u\"];\n", (void*) node, node->start, node->end);
} // draw_node


size_t
traverse(Node* const restrict node,
         size_t (*func)(Node const* const restrict, void* const restrict),
         void* const restrict arg)
{
    if (node == NULL)
    {
        return 0;
    } // if

    return traverse(node->leftsibling, func, arg) +
           func(node, arg) +
           traverse(node->rightchild, func, arg);
} // traverse


size_t
write_dot(FILE* const restrict fp, Node* const restrict root)
{
    assert(fp != NULL);

    return fputs("digraph G {\n", fp) +
           traverse(root, draw_node, fp) +
           traverse(root, draw_edge, fp) +
           fputs("}\n", fp);
} // write_dot


void
destroy(Node* const restrict node)
{
    if (node == NULL)
    {
        return;
    } // if

    destroy(node->leftsibling);
    destroy(node->rightchild);
    free(node);
} // destroy


void
read_bed(FILE* const restrict fp,
         uint32_t const max,
         Node** tree,
         Node* start_table[])
{
    assert(fp != NULL);

    // push root node on the stack
    Node* last_node = create_node(0, max);
    *tree = last_node;

    Stack* restrict stack = stack_push(NULL, last_node);
    Node* restrict top = stack_top(stack);

    // read first entry
    char reference[128] = {0};
    uint32_t start = 0;
    uint32_t end = 0;
    int ret = fscanf(fp, "%127s %" PRIu32 " %" PRIu32, reference, &start, &end);
    assert(ret == 3);

    for (uint32_t i = 0; i < max; ++i)
    {
        //fprintf(stderr, "index iteration: %" PRIu32 "\n", i);

        // pop while end smaller than index
        while (top->end < i)
        {
            assert(stack != NULL);
            //fprintf(stderr, "Popping %" PRIu32 "--%" PRIu32 " from stack\n", top->start, top->end);
            stack = stack_pop(stack);
            top = stack_top(stack);
        } // while

        while (i >= start && i <= end && !(start == top->start && end == top->end))
        {
            do
            {
                Node* const node = create_node(start, end);
                //fprintf(stderr, "Pushing %" PRIu32 "--%" PRIu32 " on stack\n", start, end);
                stack = stack_push(stack, node);
                top = stack_top(stack);

                last_node = add_node(last_node, node);

                ret = fscanf(fp, "%127s %" PRIu32 " %" PRIu32, reference, &start, &end);
                if (ret == -1)
                {
                    break;
                } // if
                assert(ret == 3);
           } while (start == top->start && end == top->end);
        } /// while

        //fprintf(stderr, "start_table[%" PRIu32 "] := %" PRIu32 "--%" PRIu32 "\n", i, top->start, top->end);
        start_table[i] = top;
    } // for

    assert(stack != NULL);
    stack = stack_pop(stack);
    assert(stack == NULL);
} // read_bed


size_t
trav(Node* const node,
     uint32_t const start,
     uint32_t const end,
     Stack* restrict* const restrict output)
{
    size_t count = 0;
    Node* restrict tmp = node;
    while (tmp != NULL && end <= tmp->end)
    {
        //fprintf(stderr, "trav: %p: %" PRIu32 "--%" PRIu32 "\n", (void*) tmp, tmp->start, tmp->end);
        if (start >= tmp->start)
        {
            count += trav(tmp->rightchild, start, end, output) + 1;
            *output = stack_push(*output, tmp);
        } // if
        tmp = tmp->leftsibling;
    } // while
    return count;
} // trav


size_t
stab(Node* const node,
     uint32_t const start,
     uint32_t const end,
     Stack* restrict* const restrict output)
{
    size_t count = 0;
    Node* restrict tmp = node;
    while (tmp != NULL && tmp->parent != NULL)
    {
        //fprintf(stderr, "stab: %p: %" PRIu32 "--%" PRIu32 "\n", (void*) tmp, tmp->start, tmp->end);
        if (start >= tmp->start)
        {
            count += 1;
            *output = stack_push(*output, tmp);
        } // if
        count += trav(tmp->leftsibling, start, end, output);
        tmp = tmp->parent;
    } // while
    return count;
} // stab


static struct Bed_entry bed[250000000];

void
query_bed(FILE* const restrict fp,
          Node* const start_table[restrict])
{
    assert(fp != NULL);

    TIC(read_bed);
    size_t idx = 0;
    while (fscanf(fp, "%31s %d %d", bed[idx].reference, &bed[idx].start, &bed[idx].end) == 3)
    {
        idx += 1;
    } // while
    TOC(read_bed);

    TIC(annotate);
    for (size_t i = 0; i < idx; ++i)
    {
        Stack* restrict output = NULL;
        bed[i].count = stab(start_table[bed[i].end], bed[i].start, bed[i].end, &output);
        bed[i].node_count = 0;
        /*
        while (output != NULL)
        {
            output = stack_pop(output);
            bed[i].count += 1;
        } // while
        */
    } // for
    TOC(annotate);

    TIC(write_ann);
    for (size_t i = 0; i < idx; ++i)
    {
        printf("%s\t%d\t%d\t%zu\t%zu\n", bed[i].reference, bed[i].start, bed[i].end, bed[i].count, bed[i].node_count);
    } // for
    TOC(write_ann);
} // query_bed