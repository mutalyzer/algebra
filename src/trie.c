#include <stddef.h>     // NULL, size_t
#include <stdlib.h>     // free, malloc
#include <string.h>     // strncpy

#include "trie.h"   // Trie, trie_*


inline Trie*
trie_init(size_t const len, char const key[static len])
{
    Trie* const trie = malloc(sizeof(*trie));
    if (trie == NULL)
    {
        return NULL;
    } // if

    trie->key = malloc(len);
    if (trie->key == NULL)
    {
        free(trie);
        return NULL;
    } // if

    strncpy(trie->key, key, len);
    trie->len = len;

    trie->link = NULL;
    trie->next = NULL;

    return trie;
} // trie_init


void
trie_destroy(Trie* const trie)
{
    if (trie == NULL)
    {
        return;
    } // if

    trie_destroy(trie->link);
    trie_destroy(trie->next);
    free(trie->key);
    free(trie);
} // if


static inline size_t
prefix_length(size_t const len_lhs, char const lhs[static len_lhs],
    size_t const len_rhs, char const rhs[static len_rhs])
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


Trie const*
trie_find(Trie* const trie, size_t const len, char const key[static len])
{
    if (trie == NULL)
    {
        return NULL;
    } // if

    size_t const k = prefix_length(len, key, trie->len, trie->key);
    if (k == 0)
    {
        return trie_find(trie->next, len, key);
    } // if
    if (k == len)
    {
        return trie;
    } // if
    if (k == trie->len)
    {
        return trie_find(trie->link, len - k, key + k);
    } // if

    return NULL;
} // trie_find


static inline void
split(Trie* const trie, size_t const k)
{
    Trie* const suffix = trie_init(trie->len - k, trie->key + k);
    if (suffix == NULL)
    {
        return;
    } // if

    char* const key = malloc(k);
    if (key == NULL)
    {
        trie_destroy(suffix);
        return;
    } // if
    strncpy(key, trie->key, k);
    free(trie->key);

    suffix->link = trie->link;
    trie->link = suffix;

    trie->key = key;
    trie->len = k;
} // split


Trie*
trie_insert(Trie* const trie, size_t const len, char const key[static len])
{
    if (trie == NULL)
    {
        return trie_init(len, key);
    } // if

    size_t const k = prefix_length(len, key, trie->len, trie->key);
    if (k == 0)
    {
        trie->next = trie_insert(trie->next, len, key);
    } // if
    else if (k < len)
    {
        if (k < trie->len)
        {
            split(trie, k);
        } // if
        trie->link = trie_insert(trie->link, len - k, key + k);
    } // if
    return trie;
} // trie_insert
