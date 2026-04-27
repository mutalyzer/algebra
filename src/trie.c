#include <stdbool.h>    // true
#include <stddef.h>     // NULL, size_t
#include <string.h>     // memcpy

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/string.h"      // GVA_String, gva_prefix_length
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "array.h"  // ARRAY_*, array_*
#include "trie.h"   // Trie, TrieNode, trie_*


inline Trie
trie_init(GVA_Allocator const allocator)
{
    return (Trie)
    {
        .allocator = allocator,
        .root = GVA_NULL,
    };
} // trie_init


inline void
trie_destroy(Trie self[static 1])
{
    self->strings = ARRAY_DESTROY(self->allocator, self->strings);
    self->nodes = ARRAY_DESTROY(self->allocator, self->nodes);
    self->root = GVA_NULL;
} // trie_destroy


static inline size_t
concat(Trie self[static restrict 1],
    size_t const len, char const key[static restrict len])
{
    size_t const start = array_length(self->strings);
    self->strings = array_ensure(self->allocator, self->strings, sizeof(*self->strings), len);
    if (self->strings == NULL)
    {
        return start;  // OOM
    } // if
    memcpy(self->strings + start, key, len);
    array_header(self->strings)->length += len;
    return start;
} // concat


gva_uint
trie_insert(Trie self[static restrict 1],
    size_t const len, char const key[static restrict len])
{
    if (self->root == GVA_NULL)
    {
        gva_uint const start = concat(self, len, key);
        self->root = ARRAY_APPEND(self->allocator, self->nodes,
            ((TrieNode)
            {
                .link = GVA_NULL,
                .next = GVA_NULL,
                .p_start = start,
                .start = start,
                .end = len,
            }));
        return self->root;
    } // if

    gva_uint prev = GVA_NULL;
    size_t prefix = 0;
    gva_uint idx = self->root;
    while (true)
    {
        gva_uint const p_len = self->nodes[idx].end - self->nodes[idx].p_start;
        size_t const k = gva_prefix_length(len - prefix, key + prefix,
            p_len, self->strings + self->nodes[idx].p_start);
        if (k == len - prefix && k == p_len)
        {
            return idx;  // existing found
        } // if

        if (k == 0)
        {
            if (self->nodes[idx].next == GVA_NULL)
            {
                gva_uint const start = concat(self, len, key);
                gva_uint const next = ARRAY_APPEND(self->allocator, self->nodes,
                    ((TrieNode)
                    {
                        .link = GVA_NULL,
                        .next = GVA_NULL,
                        .p_start = start + prefix,
                        .start = start,
                        .end = array_length(self->strings),
                    }));
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
                gva_uint const start = concat(self, len, key);
                gva_uint const link = ARRAY_APPEND(self->allocator, self->nodes,
                    ((TrieNode)
                    {
                        .link = GVA_NULL,
                        .next = GVA_NULL,
                        .p_start = start + prefix,
                        .start = start,
                        .end = array_length(self->strings),
                    }));
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
                gva_uint const start = concat(self, len, key);
                next = ARRAY_APPEND(self->allocator, self->nodes,
                    ((TrieNode)
                    {
                        .link = GVA_NULL,
                        .next = GVA_NULL,
                        .p_start = start + prefix,
                        .start = start,
                        .end = array_length(self->strings),
                    }));
            } // if
            gva_uint const link = ARRAY_APPEND(self->allocator, self->nodes,
                ((TrieNode)
                {
                    .link = idx,
                    .next = self->nodes[idx].next,
                    .p_start = self->nodes[idx].p_start,
                    .start = self->nodes[idx].start,
                    .end = self->nodes[idx].p_start + k,
                }));
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


inline GVA_String
trie_string(Trie const self, size_t const idx)
{
    if (idx >= array_length(self.nodes))
    {
        return (GVA_String) {0};
    } // if
    return (GVA_String)
    {
        .len = self.nodes[idx].end - self.nodes[idx].start,
        .str = self.strings + self.nodes[idx].start,
    };
} //trie_string
