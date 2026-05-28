#ifndef GVA_TRIE_H
#define GVA_TRIE_H


#include <stddef.h>     // size_t

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/string.h"      // GVA_String
#include "../include/types.h"       // gva_uint


typedef struct
{
    gva_uint link;
    gva_uint next;
    gva_uint p_start;
    gva_uint start;
    gva_uint end;
} Trie_Node;


typedef struct
{
    char* restrict      strings;
    Trie_Node* restrict nodes;
    GVA_Allocator       nodes_allocator;
    GVA_Allocator       strings_allocator;
    gva_uint            root;
} Trie;


Trie
trie_init(GVA_Allocator const nodes_allocator, GVA_Allocator const strings_allocator);


void
trie_destroy(Trie self[static 1]);


gva_uint
trie_insert(Trie self[static restrict 1],
    size_t const len, char const key[static restrict len]);


GVA_String
trie_string(Trie const self, size_t const idx);


#endif  // GVA_TRIE_H
