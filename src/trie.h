#ifndef GVA_TRIE_H
#define GVA_TRIE_H


#include <stddef.h>     // size_t

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/types.h"       // gva_uint


typedef struct
{
    gva_uint link;
    gva_uint next;
    gva_uint p_start;
    gva_uint start;
    gva_uint end;
} TrieNode;


typedef struct
{
    char*     strings;
    TrieNode* nodes;
    gva_uint  root;
} Trie;


gva_uint
trie_insert(GVA_Allocator const allocator, Trie self[static restrict 1],
    size_t const len, char const key[static restrict len]);


void
trie_destroy(GVA_Allocator const allocator, Trie self[static 1]);


#endif  // GVA_TRIE_H
