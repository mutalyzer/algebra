#include <stdbool.h>    // true
#include <stddef.h>     // size_t
#include <string.h>     // memcpy

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "array.h"  // ARRAY_*, array_*
#include "trie.h"   // Trie, TrieNode, trie_*


inline void
trie_destroy(GVA_Allocator const allocator, Trie self[static 1])
{
    self->strings = ARRAY_DESTROY(allocator, self->strings);
    self->nodes = ARRAY_DESTROY(allocator, self->nodes);
} // trie_destroy


static inline size_t
prefix_length(size_t const len_lhs, char const lhs[static restrict len_lhs],
    size_t const len_rhs, char const rhs[static restrict len_rhs])
{
    for (size_t i = 0; i < len_lhs; ++i)
    {
        if (i >= len_rhs || lhs[i] != rhs[i])
        {
            return i;
        } // if
    } // for
    return len_lhs;
} // prefix_length


static inline size_t
concat(GVA_Allocator const allocator, Trie self[static restrict 1],
    size_t const len, char const key[static restrict len])
{
    size_t const start = array_length(self->strings);
    self->strings = array_ensure(allocator, self->strings, sizeof(*self->strings), len);
    memcpy(self->strings + start, key, len);
    array_header(self->strings)->length += len;
    return start;
} // concat


gva_uint
trie_insert(GVA_Allocator const allocator, Trie self[static restrict 1],
    size_t const len, char const key[static restrict len])
{
    if (self->root == GVA_NULL)
    {
        gva_uint const start = concat(allocator, self, len, key);
        self->root = ARRAY_APPEND(allocator, self->nodes, ((TrieNode) {GVA_NULL, GVA_NULL, start, start, len})) - 1;
        return self->root;
    } // if

    gva_uint prev = GVA_NULL;
    size_t prefix = 0;
    gva_uint idx = self->root;
    while (true)
    {
        gva_uint const p_len = self->nodes[idx].end - self->nodes[idx].p_start;
        size_t const k = prefix_length(len - prefix, key + prefix, p_len, self->strings + self->nodes[idx].p_start);
        if (k == len - prefix && k == p_len)
        {
            return idx;  // existing found
        } // if

        if (k == 0)
        {
            if (self->nodes[idx].next == GVA_NULL)
            {
                gva_uint const start = concat(allocator, self, len, key);
                gva_uint const end = array_length(self->strings);
                gva_uint const next = ARRAY_APPEND(allocator, self->nodes, ((TrieNode) {GVA_NULL, GVA_NULL, start + prefix, start, end})) - 1;
                self->nodes[idx].next = next;
                return next;
            } // if
            prev = idx;
            idx = self->nodes[idx].next;
        } // if
        else if (k == p_len)
        {
            prefix += k;
            if (self->nodes[idx].link == GVA_NULL)
            {
                gva_uint const start = concat(allocator, self, len, key);
                gva_uint const end = array_length(self->strings);
                gva_uint const link = ARRAY_APPEND(allocator, self->nodes, ((TrieNode) {GVA_NULL, GVA_NULL, start + prefix, start, end})) - 1;
                self->nodes[idx].link = link;
                return link;
            } // if
            prev = idx;
            idx = self->nodes[idx].link;
        } // if
        else
        {
            prefix += k;
            gva_uint next = GVA_NULL;
            if (prefix < len)
            {
                gva_uint const start = concat(allocator, self, len, key);
                gva_uint const end = array_length(self->strings);
                next = ARRAY_APPEND(allocator, self->nodes, ((TrieNode) {GVA_NULL, GVA_NULL, start + prefix, start, end})) - 1;
            } // if
            gva_uint const link = ARRAY_APPEND(allocator, self->nodes, ((TrieNode) {idx, self->nodes[idx].next, self->nodes[idx].p_start, self->nodes[idx].start, self->nodes[idx].p_start + k})) - 1;
            self->nodes[idx].next = next;
            self->nodes[idx].p_start += k;
            if (prev != GVA_NULL)
            {
                if (self->nodes[prev].link == idx)
                {
                    self->nodes[prev].link = link;
                } // if
                else
                {
                    self->nodes[prev].next = link;
                } // else
            } // if
            else
            {
                self->root = link;
            } // else
            if (next != GVA_NULL)
            {
                return next;
            } // if
            return link;
        } // else
    } // while
} // trie_insert