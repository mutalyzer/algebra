#ifndef GVA_TRIE_H
#define GVA_TRIE_H


#include <stddef.h>     // size_t


typedef struct Trie
{
    struct Trie* link;
    struct Trie* next;
    size_t len;
    char*  key;
} Trie;


Trie*
trie_init(size_t const len, char const key[static len]);


void
trie_destroy(Trie* const trie);


Trie const*
trie_find(Trie* const trie, size_t const len, char const key[static len]);


Trie*
trie_insert(Trie* const trie, size_t const len, char const key[static len]);


#endif  // GVA_TRIE_H
