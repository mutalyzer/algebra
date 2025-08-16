#include <stddef.h>     // NULL
#include <stdint.h>     // uint64_t

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/types.h"       // GVA_NULL
#include "array.h"          // ARRAY_DESTROY
#include "common.h"         // MAX
#include "interval_tree.h"  // Interval_Tree, interval_tree_*


inline Interval_Tree
interval_tree_init(void)
{
    return (Interval_Tree) {NULL, GVA_NULL};
} // interval_tree_init


inline void
interval_tree_destroy(GVA_Allocator const allocator, Interval_Tree self[static 1])
{
    self->nodes = ARRAY_DESTROY(allocator, self->nodes);
    self->root = GVA_NULL;
} // interval_tree_destroy


enum
{
    LEFT  = 0,
    RIGHT = 1
};


static inline gva_uint
calculate_max(Interval_Tree_Node const nodes[static 1], gva_uint const idx)
{
    gva_uint result = nodes[idx].end;
    if (nodes[idx].child[LEFT] != GVA_NULL)
    {
        result = MAX(result, nodes[nodes[idx].child[LEFT]].max);
    } // if
    if (nodes[idx].child[RIGHT] != GVA_NULL)
    {
        result = MAX(result, nodes[nodes[idx].child[RIGHT]].max);
    } // if
    return result;
} // calculate_max


gva_uint
interval_tree_insert(Interval_Tree self[static 1], gva_uint const idx)
{
    if (self->root == GVA_NULL)
    {
        self->root = idx;
        return idx;
    } // if

    // limiting to height 64 allows for a minimum of
    // 27,777,890,035,287 nodes
    uint64_t path = 0;  // bit-path to first unbalanced ancestor
    int len = 0;        // length of the path
    int dir = LEFT;

    gva_uint tmp = self->root;
    gva_uint tmp_par = self->root;  // parent of tmp

    gva_uint unbal = self->root;      // first unbalanced ancestor of tmp
    gva_uint unbal_par = self->root;  // parent of unbalanced

    // Insert a new node at the BST position
    while (tmp != GVA_NULL)
    {
        self->nodes[tmp].max = MAX(self->nodes[tmp].max, self->nodes[idx].end);

        if (self->nodes[tmp].balance != 0)
        {
            // this is now the first unbalanced ancestor of tmp
            unbal_par = tmp_par;
            unbal = tmp;
            path = 0;
            len = 0;
        } // if

        dir = self->nodes[idx].start > self->nodes[tmp].start;
        if (dir == RIGHT)
        {
            path |= (uint64_t) RIGHT << len;
        } // if
        len += 1;

        if (self->nodes[tmp_par].start == self->nodes[idx].start &&
            self->nodes[tmp_par].end == self->nodes[idx].end &&
            self->nodes[tmp_par].inserted == self->nodes[idx].inserted)
        {
            return tmp_par;
        } // if

        tmp_par = tmp;
        tmp = self->nodes[tmp].child[dir];
    } // while


    self->nodes[tmp_par].child[dir] = idx;

    // Update the balance factors along the path from the first
    // unbalanced ancestor to the new node
    tmp = unbal;
    while (tmp != idx)
    {
        if ((path & RIGHT) == RIGHT)
        {
            self->nodes[tmp].balance += 1;
        } // if
        else
        {
            self->nodes[tmp].balance -= 1;
        } // else

        tmp = self->nodes[tmp].child[path & RIGHT];
        path >>= 1;
    } // while

    // Do the rotations if necessary
    gva_uint root = 0;
    if (self->nodes[unbal].balance == -2)
    {
        gva_uint const child = self->nodes[unbal].child[LEFT];
        if (self->nodes[child].balance == -1)
        {
            root = child;
            self->nodes[unbal].child[LEFT] = self->nodes[child].child[RIGHT];
            self->nodes[child].child[RIGHT] = unbal;
            self->nodes[child].balance = 0;
            self->nodes[unbal].balance = 0;
            self->nodes[child].max = calculate_max(self->nodes, child);
            self->nodes[unbal].max = calculate_max(self->nodes, unbal);
        } // if
        else
        {
            root = self->nodes[child].child[RIGHT];
            self->nodes[child].child[RIGHT] = self->nodes[root].child[LEFT];
            self->nodes[root].child[LEFT] = child;
            self->nodes[unbal].child[LEFT] = self->nodes[root].child[RIGHT];
            self->nodes[root].child[RIGHT] = unbal;
            if (self->nodes[root].balance == -1)
            {
                self->nodes[child].balance = 0;
                self->nodes[unbal].balance = 1;
            } // if
            else if (self->nodes[root].balance == 0)
            {
                self->nodes[child].balance = 0;
                self->nodes[unbal].balance = 0;
            } // if
            else
            {
                self->nodes[child].balance = -1;
                self->nodes[unbal].balance = 0;
            } // else
            self->nodes[root].balance = 0;

            self->nodes[root].max = calculate_max(self->nodes, root);
            self->nodes[child].max = calculate_max(self->nodes, child);
            self->nodes[unbal].max = calculate_max(self->nodes, unbal);
        } // else
    } // if
    else if (self->nodes[unbal].balance == 2)
    {
        gva_uint const child = self->nodes[unbal].child[RIGHT];
        if (self->nodes[child].balance == 1)
        {
            root = child;
            self->nodes[unbal].child[RIGHT] = self->nodes[child].child[LEFT];
            self->nodes[child].child[LEFT] = unbal;
            self->nodes[child].balance = 0;
            self->nodes[unbal].balance = 0;

            self->nodes[child].max = calculate_max(self->nodes, child);
            self->nodes[unbal].max = calculate_max(self->nodes, unbal);
        } // if
        else
        {
            root = self->nodes[child].child[LEFT];
            self->nodes[child].child[LEFT] = self->nodes[root].child[RIGHT];
            self->nodes[root].child[RIGHT] = child;
            self->nodes[unbal].child[RIGHT] = self->nodes[root].child[LEFT];
            self->nodes[root].child[LEFT] = unbal;
            if (self->nodes[root].balance == 1)
            {
                self->nodes[child].balance = 0;
                self->nodes[unbal].balance = -1;
            } // if
            else if(self->nodes[root].balance == 0)
            {
                self->nodes[child].balance = 0;
                self->nodes[unbal].balance = 0;
            } // if
            else
            {
                self->nodes[child].balance = 1;
                self->nodes[unbal].balance = 0;
            } // else
            self->nodes[root].balance = 0;

            self->nodes[root].max = calculate_max(self->nodes, root);
            self->nodes[child].max = calculate_max(self->nodes, child);
            self->nodes[unbal].max = calculate_max(self->nodes, unbal);
        } // else
    } // if
    else
    {
        return idx;
    } // else

    if (self->root == unbal)
    {
        self->root = root;
        return idx;
    } // if

    self->nodes[unbal_par].child[unbal != self->nodes[unbal_par].child[LEFT]] = root;
    return idx;
} // interval_tree_insert
