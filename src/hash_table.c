#include <stdbool.h>    // true
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint64_t, uintmax_t, UINT64_C
#include <string.h>     // memcpy

#include "array.h"          // ARRAY_DESTROY, Array, array_*
#include "hash_table.h"     // HASH_TABLE_NOT_FOUND, hash_table_*


void*
hash_table_init(GVA_Allocator const allocator, size_t const capacity, size_t const item_size)
{
    void* const hash_table = array_init(allocator, capacity, item_size);
    if (hash_table == NULL)
    {
        return NULL;  // OOM
    } // if

    for (size_t i = 0; i < capacity; ++i)
    {
        *(uint64_t*) ((uintmax_t) hash_table + i * item_size) = HASH_TABLE_NOT_FOUND;
    } // if
    return hash_table;
} // hash_table_init


// murmur
static inline uint64_t
hash(uint64_t key)
{
    key ^= key >> 33;
    key *= UINT64_C(0xff51afd7ed558ccd);
    key ^= key >> 33;
    key *= UINT64_C(0xc4ceb9fe1a85ec53);
    key ^= key >> 33;
    return key;
} // hash


size_t
hash_table_index(void* const self, size_t const item_size, size_t const key)
{
    size_t probe = 0;
    while (true)
    {
        // Quadratic probing with c1 and c2 1/2.
        size_t const idx = (hash(key) + (probe + probe * probe) / 2) % array_header(self)->capacity;
        if (*(uint64_t*) ((uintmax_t) self + idx * item_size) == HASH_TABLE_NOT_FOUND ||
            *(uint64_t*) ((uintmax_t) self + idx * item_size) == key)
        {
            return idx;
        } // if
        probe += 1;
    } // while
} // hash_table_index


void*
hash_table_ensure(GVA_Allocator const allocator, void* const self, size_t const item_size)
{
    size_t const old_capacity = array_header(self)->capacity;
    if (array_length(self) < old_capacity / 2)
    {
        return self;
    } // if
    void* new_table = hash_table_init(allocator, old_capacity * 2, item_size);  // OVERFLOW
    if (new_table == NULL)
    {
        return self;  // OOM
    } // if
    for (size_t i = 0; i < old_capacity; ++i)
    {
        void const* const src = (void*) ((uintmax_t) self + i * item_size);
        if (*(uint64_t*) src != HASH_TABLE_NOT_FOUND)
        {
            array_header(new_table)->length += 1;
            memcpy((void*) ((uintmax_t) new_table + hash_table_index(new_table, item_size, *(uint64_t*) src) * item_size),
                src, item_size);
        } // if
    } // for
    ARRAY_DESTROY(allocator, (Array*) self);
    return new_table;
} // hash_table_ensure
