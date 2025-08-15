#include <stdint.h>     // uint64_t
#include <stddef.h>     // NULL, size_t

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/index.h"       // GVA_Index
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "array.h"      // ARRAY_APPEND, ARRAY_DESTROY, array_length
#include "common.h"     // MAX


typedef struct Allele
{
    size_t   data;
    gva_uint start;
    gva_uint end;
} Allele;


typedef struct Node
{
    gva_uint child[2];
    gva_uint start;
    gva_uint end;
    gva_uint max;
    gva_uint allele;
    // FIXME: we assume gva_uint is 32bits
    gva_uint balance   :  3;
    gva_uint inserted  : 29;
} Node;


typedef struct Next_Allele
{
    gva_uint allele;
    gva_uint next;
} Next_Allele;


inline GVA_Index
gva_index_init(GVA_Allocator const allocator)
{
    return (GVA_Index) {allocator, GVA_NULL, NULL, NULL, NULL};
} // gva_index_init


inline void
gva_index_destroy(GVA_Index self)
{
    self.root = GVA_NULL;
    self.alleles = ARRAY_DESTROY(self.allocator, self.alleles);
    self.nodes = ARRAY_DESTROY(self.allocator, self.nodes);
    self.next_alleles = ARRAY_DESTROY(self.allocator, self.next_alleles);
} // gva_index_destroy


enum
{
    LEFT  = 0,
    RIGHT = 1
};


static inline gva_uint
update_max(GVA_Index const self, gva_uint const idx)
{
    gva_uint result = self.nodes[idx].end;
    if (self.nodes[idx].child[LEFT] != GVA_NULL)
    {
        result = MAX(result, self.nodes[self.nodes[idx].child[LEFT]].max);
    } // if
    if (self.nodes[idx].child[RIGHT] != GVA_NULL)
    {
        result = MAX(result, self.nodes[self.nodes[idx].child[RIGHT]].max);
    } // if
    return result;
} // update_max


static gva_uint
insert(GVA_Index self, gva_uint const idx)
{
    if (self.root == GVA_NULL)
    {
        self.root = idx;
        return idx;
    } // if

    // limiting to height 64 becomes a problem after allocating 413 TiB
    // at the earliest; it allows for a minimum of
    // 27,777,890,035,287 nodes
    uint64_t path = 0;  // bit-path to first unbalanced ancestor
    int len = 0;        // length of the path
    int dir = 0;

    gva_uint tmp = self.root;
    gva_uint tmp_par = self.root;  // parent of tmp

    gva_uint unbal = self.root;      // first unbalanced ancestor of tmp
    gva_uint unbal_par = self.root;  // parent of unbalanced

    // Insert a new node at the BST position
    while (tmp != GVA_NULL)
    {
        self.nodes[tmp].max = MAX(self.nodes[tmp].max, self.nodes[idx].end);

        if (self.nodes[tmp].balance != 0)
        {
            // this is now the first unbalanced ancestor of tmp
            unbal_par = tmp_par;
            unbal = tmp;
            path = 0;
            len = 0;
        } // if

        dir = self.nodes[idx].start > self.nodes[tmp].start;
        if (dir == RIGHT)
        {
            path |= (uint64_t) RIGHT << len;
        } // if
        len += 1;

        tmp_par = tmp;
        tmp = self.nodes[tmp].child[dir];
    } // while

    // TODO: early return if found

    self.nodes[tmp_par].child[dir] = idx;

    // Update the balance factors along the path from the first
    // unbalanced ancestor to the new node
    tmp = unbal;
    while (tmp != idx)
    {
        if ((path & RIGHT) == RIGHT)
        {
            self.nodes[tmp].balance += 1;
        } // if
        else
        {
            self.nodes[tmp].balance -= 1;
        } // else

        tmp = self.nodes[tmp].child[path & RIGHT];
        path >>= 1;
    } // while

    // Do the rotations if necessary
    gva_uint root = 0;
    if (self.nodes[unbal].balance == -2)
    {
        gva_uint const child = self.nodes[unbal].child[LEFT];
        if (self.nodes[child].balance == -1)
        {
            root = child;
            self.nodes[unbal].child[LEFT] = self.nodes[child].child[RIGHT];
            self.nodes[child].child[RIGHT] = unbal;
            self.nodes[child].balance = 0;
            self.nodes[unbal].balance = 0;
            self.nodes[child].max = update_max(self, child);
            self.nodes[unbal].max = update_max(self, unbal);
        } // if
        else
        {
            root = self.nodes[child].child[RIGHT];
            self.nodes[child].child[RIGHT] = self.nodes[root].child[LEFT];
            self.nodes[root].child[LEFT] = child;
            self.nodes[unbal].child[LEFT] = self.nodes[root].child[RIGHT];
            self.nodes[root].child[RIGHT] = unbal;
            if (self.nodes[root].balance == -1)
            {
                self.nodes[child].balance = 0;
                self.nodes[unbal].balance = 1;
            } // if
            else if (self.nodes[root].balance == 0)
            {
                self.nodes[child].balance = 0;
                self.nodes[unbal].balance = 0;
            } // if
            else
            {
                self.nodes[child].balance = -1;
                self.nodes[unbal].balance = 0;
            } // else
            self.nodes[root].balance = 0;

            self.nodes[root].max = update_max(self, root);
            self.nodes[child].max = update_max(self, child);
            self.nodes[unbal].max = update_max(self, unbal);
        } // else
    } // if
    else if (self.nodes[unbal].balance == 2)
    {
        gva_uint const child = self.nodes[unbal].child[RIGHT];
        if (self.nodes[child].balance == 1)
        {
            root = child;
            self.nodes[unbal].child[RIGHT] = self.nodes[child].child[LEFT];
            self.nodes[child].child[LEFT] = unbal;
            self.nodes[child].balance = 0;
            self.nodes[unbal].balance = 0;

            self.nodes[child].max = update_max(self, child);
            self.nodes[unbal].max = update_max(self, unbal);
        } // if
        else
        {
            root = self.nodes[child].child[LEFT];
            self.nodes[child].child[LEFT] = self.nodes[root].child[RIGHT];
            self.nodes[root].child[RIGHT] = child;
            self.nodes[unbal].child[RIGHT] = self.nodes[root].child[LEFT];
            self.nodes[root].child[LEFT] = unbal;
            if (self.nodes[root].balance == 1)
            {
                self.nodes[child].balance = 0;
                self.nodes[unbal].balance = -1;
            } // if
            else if(self.nodes[root].balance == 0)
            {
                self.nodes[child].balance = 0;
                self.nodes[unbal].balance = 0;
            } // if
            else
            {
                self.nodes[child].balance = 1;
                self.nodes[unbal].balance = 0;
            } // else
            self.nodes[root].balance = 0;

            self.nodes[root].max = update_max(self, root);
            self.nodes[child].max = update_max(self, child);
            self.nodes[unbal].max = update_max(self, unbal);
        } // else
    } // if
    else
    {
        return idx;
    } // else

    if (self.root == unbal)
    {
        self.root = root;
        return idx;
    } // if

    self.nodes[unbal_par].child[unbal != self.nodes[unbal_par].child[LEFT]] = root;
    return idx;
} // insert


void
gva_index_add(GVA_Index self, size_t const data,
    size_t const n, GVA_Variant const variants[static n])
{
    gva_uint const allele = ARRAY_APPEND(self.allocator, self.alleles, ((Allele) {data, array_length(self.next_alleles) - 1, array_length(self.next_alleles) + n - 1})) - 1;
    for (size_t i = 0; i < n; ++i)
    {
        gva_uint const idx = ARRAY_APPEND(self.allocator, self.next_alleles, ((Next_Allele) {allele, GVA_NULL})) - 1;
        if (1 /* not found */)
        {
            gva_uint const inserted = 0; /* insert in trie */
            gva_uint const node = ARRAY_APPEND(self.allocator, self.nodes, ((Node) {{GVA_NULL, GVA_NULL}, variants[i].start, variants[i].end, variants[i].end, idx, 0, inserted})) - 1;
            insert(self, node);
        } // if
        else
        {
            self.next_alleles[idx].next = 0; /* found node */
            self.nodes[0 /* found node */].allele = idx;
        } // else
    } // for
} // gva_index_add
