#include <stddef.h>     // NULL
#include <stdint.h>     // uint8_t

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
update_max(Interval_Tree const self[static 1], gva_uint const idx)
{
    gva_uint result = self->nodes[idx].end;
    if (self->nodes[idx].child[LEFT] != GVA_NULL)
    {
        result = MAX(result, self->nodes[self->nodes[idx].child[LEFT]].max);
    } // if
    if (self->nodes[idx].child[RIGHT] != GVA_NULL)
    {
        result = MAX(result, self->nodes[self->nodes[idx].child[RIGHT]].max);
    } // if
    return result;
} // update_max


static int
interval_cmp(gva_uint const lhs_start, gva_uint const lhs_end, gva_uint const lhs_inserted,
    gva_uint const rhs_start, gva_uint const rhs_end, gva_uint const rhs_inserted)
{
    if (lhs_start < rhs_start)
    {
        return LEFT;
    } // if
    if (rhs_start < lhs_start)
    {
        return RIGHT;
    } // if
    if (lhs_end < rhs_end)
    {
        return LEFT;
    } // if
    if (rhs_end < lhs_end)
    {
        return RIGHT;
    } // if
    if (lhs_inserted < rhs_inserted)
    {
        return LEFT;
    } // if
    if (rhs_inserted < lhs_inserted)
    {
        return RIGHT;
    } // if
    return -1;  // identity
} // interval_cmp


gva_uint
interval_tree_insert(Interval_Tree self[static 1], gva_uint const idx)
{
    if (self->root == GVA_NULL)
    {
        self->root = idx;
        return idx;
    } // if

    gva_uint y = self->root;  // Top node to update balance factor,
    gva_uint z = self->root;  // and parent.

    uint8_t path[64];  // Cached comparison results.
    int k = 0;         // Number of cached results.
    int dir = LEFT;    // Direction to descend.

    gva_uint p = y;  // Iterator,
    gva_uint q = z;  // and parent.
    while (p != GVA_NULL)
    {
        dir = interval_cmp(self->nodes[idx].start, self->nodes[idx].end, self->nodes[idx].inserted,
            self->nodes[p].start, self->nodes[p].end, self->nodes[p].inserted);
        if (dir == -1)
        {
            return p;  // found
        } // if

        self->nodes[p].max = MAX(self->nodes[p].max, self->nodes[idx].end);

        if (self->nodes[p].balance != 0)
        {
            z = q;
            y = p;
            k = 0;
        } // if
        path[k] = dir;
        k += 1;
        q = p;
        p = self->nodes[p].child[dir];
    } // while

    // add
    self->nodes[q].child[dir] = idx;

    // update balance factors
    p = y;
    k = 0;
    while (p != idx)
    {
        if (path[k] == LEFT)
        {
            self->nodes[p].balance -= 1;
        } // if
        else
        {
            self->nodes[p].balance += 1;
        } // else
        p = self->nodes[p].child[path[k]];
        k += 1;
    } // for

    // rotations
    gva_uint w;  // New root of rebalanced subtree.
    if (self->nodes[y].balance == -2)
    {
        gva_uint const x = self->nodes[y].child[LEFT];
        if (self->nodes[x].balance == -1)
        {
            w = x;
            self->nodes[y].child[LEFT] = self->nodes[x].child[RIGHT];
            self->nodes[x].child[RIGHT] = y;
            self->nodes[x].balance = 0;
            self->nodes[y].balance = 0;

            self->nodes[x].max = update_max(self, x);
            self->nodes[y].max = update_max(self, y);
        } // if
        else
        {
            w = self->nodes[x].child[RIGHT];
            self->nodes[x].child[RIGHT] = self->nodes[w].child[LEFT];
            self->nodes[w].child[LEFT] = x;
            self->nodes[y].child[LEFT] = self->nodes[w].child[RIGHT];
            self->nodes[w].child[RIGHT] = y;
            if (self->nodes[w].balance == -1)
            {
                self->nodes[x].balance = 0;
                self->nodes[y].balance = +1;
            } // if
            else if (self->nodes[w].balance == 0)
            {
                self->nodes[x].balance = 0;
                self->nodes[y].balance = 0;
            } // if
            else
            {
                self->nodes[x].balance = -1;
                self->nodes[y].balance = 0;
            } // else
            self->nodes[w].balance = 0;

            self->nodes[x].max = update_max(self, x);
            self->nodes[y].max = update_max(self, y);
            self->nodes[w].max = update_max(self, w);
        } // else
    } // if
    else if (self->nodes[y].balance == +2)
    {
        gva_uint const x = self->nodes[y].child[RIGHT];
        if (self->nodes[x].balance == +1)
        {
            w = x;
            self->nodes[y].child[RIGHT] = self->nodes[x].child[LEFT];
            self->nodes[x].child[LEFT] = y;
            self->nodes[x].balance = 0;
            self->nodes[y].balance = 0;

            self->nodes[x].max = update_max(self, x);
            self->nodes[y].max = update_max(self, y);
        } // if
        else
        {
            w = self->nodes[x].child[LEFT];
            self->nodes[x].child[LEFT] = self->nodes[w].child[RIGHT];
            self->nodes[w].child[RIGHT] = x;
            self->nodes[y].child[RIGHT] = self->nodes[w].child[LEFT];
            self->nodes[w].child[LEFT] = y;
            if (self->nodes[w].balance == +1)
            {
                self->nodes[x].balance = 0;
                self->nodes[y].balance = -1;
            } // if
            else if (self->nodes[w].balance == 0)
            {
                self->nodes[x].balance = 0;
                self->nodes[y].balance = 0;
            } // if
            else
            {
                self->nodes[x].balance = +1;
                self->nodes[y].balance = 0;
            } // else
            self->nodes[w].balance = 0;

            self->nodes[x].max = update_max(self, x);
            self->nodes[y].max = update_max(self, y);
            self->nodes[w].max = update_max(self, w);
        } // else
    } // if
    else
    {
        return idx;
    } // else

    if (self->root == y)
    {
        self->root = w;
        return idx;
    } // if

    self->nodes[z].child[y != self->nodes[z].child[LEFT]] = w;
    return idx;
} // interval_tree_insert


// FIXME: recursion!
static gva_uint*
intersect(GVA_Allocator const allocator, Interval_Tree_Node const nodes[static 1],
    gva_uint const idx, gva_uint const start, gva_uint const end, gva_uint* results)
{
    if (idx == GVA_NULL || start > nodes[idx].max)
    {
        return results;
    } // if

    results = intersect(allocator, nodes, nodes[idx].child[LEFT], start, end, results);

    if (end < nodes[idx].start)
    {
        return results;
    } // if

    if (start <= nodes[idx].end)  // already satisfied above: && end >= nodes[idx].start)
    {
        ARRAY_APPEND(allocator, results, idx);
    } // if

    return intersect(allocator, nodes, nodes[idx].child[RIGHT], start, end, results);
} // intersect


inline gva_uint*
interval_tree_intersection(GVA_Allocator const allocator, Interval_Tree const self,
    gva_uint const start, gva_uint const end)
{
    return intersect(allocator, self.nodes, self.root, start, end, NULL);
} // interval_tree_intersection
