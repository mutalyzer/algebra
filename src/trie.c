#include <stdbool.h>    // true
#include <stddef.h>     // size_t
#include <string.h>     // memcpy

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/types.h"       // GVA_NULL
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


void
trie_insert(GVA_Allocator const allocator, Trie self[static restrict 1],
    size_t const len, char const key[static restrict len])
{
    if (self->nodes == NULL)
    {
        ARRAY_APPEND(allocator, self->nodes, ((TrieNode) {GVA_NULL, GVA_NULL, concat(allocator, self, len, key), len}));
        return;
    } // if

    size_t prefix = 0;
    gva_uint idx = 0;
    while (true)
    {
        size_t const k = prefix_length(len - prefix, key + prefix, self->nodes[idx].end - self->nodes[idx].start, self->strings + self->nodes[idx].start);
        if (k == len - prefix && k == self->nodes[idx].end - self->nodes[idx].start)
        {
            return;  // duplicate found
        } // if

        if (k == 0)
        {
            if (self->nodes[idx].next == GVA_NULL)
            {
                gva_uint const next = ARRAY_APPEND(allocator, self->nodes, ((TrieNode) {GVA_NULL, GVA_NULL, concat(allocator, self, len - prefix, key + prefix), array_length(self->strings)})) - 1;
                self->nodes[idx].next = next;
                return;
            } // if
            idx = self->nodes[idx].next;
        } // if
        else if (k == self->nodes[idx].end - self->nodes[idx].start)
        {
            prefix += k;
            if (self->nodes[idx].link == GVA_NULL)
            {
                gva_uint const link = ARRAY_APPEND(allocator, self->nodes, ((TrieNode) {GVA_NULL, GVA_NULL, concat(allocator, self, len - prefix, key + prefix), array_length(self->strings)})) - 1;
                self->nodes[idx].link = link;
                return;
            } // if
            idx = self->nodes[idx].link;
        } // if
        else
        {
            gva_uint next = GVA_NULL;
            if (k < len - prefix)
            {
                prefix += k;
                next = ARRAY_APPEND(allocator, self->nodes, ((TrieNode) {GVA_NULL, GVA_NULL, concat(allocator, self, len - prefix, key + prefix), array_length(self->strings)})) - 1;
            } // if
            gva_uint const link = ARRAY_APPEND(allocator, self->nodes, ((TrieNode) {self->nodes[idx].link, next, self->nodes[idx].start + k, self->nodes[idx].end})) - 1;
            self->nodes[idx].link = link;
            self->nodes[idx].end = self->nodes[idx].start + k;
            return;
        } // else
    } // while
} // trie_insert

