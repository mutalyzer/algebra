#include <stdbool.h>    // true
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint32_t, uintmax_t
#include <string.h>     // memcpy

#include "array.h"          // ARRAY_DESTROY, Array, array_*
#include "hash_table.h"     // hash_table_*


void*
hash_table_init(GVA_Allocator const allocator, size_t const capacity, size_t const item_size)
{
    void* hash_table = array_init(allocator, capacity, item_size);
    for (size_t i = 0; i < capacity; ++i)
    {
        *(uint32_t*) ((uintmax_t) hash_table + i * item_size) = -1;
    } // if
    return hash_table;
} // hash_table_init


// Jenkins hash
static inline uint32_t
hash(uint32_t key)
{
    key = (key + 0x7ed55d16) + (key << 12);
    key = (key ^ 0xc761c23c) ^ (key >> 19);
    key = (key + 0x165667b1) + (key << 5);
    key = (key + 0xd3a2646c) ^ (key << 9);
    key = (key + 0xfd7046c5) + (key << 3);
    key = (key ^ 0xb55a4f09) ^ (key >> 16);
    return key;
} // hash


size_t
hash_table_index(void* const self, size_t const item_size, uint32_t const key)
{
    size_t probe = 0;
    while (true)
    {
        // Quadratic probing with c1 and c2 1/2.
        size_t const idx = (hash(key) + (probe + probe * probe) / 2) % array_header(self)->capacity;
        if (*(uint32_t*) ((uintmax_t) self + idx * item_size) == (uint32_t) -1 ||
            *(uint32_t*) ((uintmax_t) self + idx * item_size) == key)
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
        if (*(uint32_t*) src != (uint32_t) -1)
        {
            array_header(new_table)->length += 1;
            memcpy(
                (void*) ((uintmax_t) new_table + hash_table_index(new_table, item_size, *(uint32_t*) src) * item_size),
                src, item_size);
        } // if
    } // for
    ARRAY_DESTROY(allocator, (Array*) self);
    return new_table;
} // hash_table_ensure
