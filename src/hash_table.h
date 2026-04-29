#ifndef GVA_HASH_TABLE_H
#define GVA_HASH_TABLE_H


#include <stddef.h>     // size_t
#include <stdint.h>     // uint64_t, UINT64_MAX

#include "../include/allocator.h"   // GVA_Allocator
#include "array.h"      // ARRAY_DESTROY, array_header


#define HASH_TABLE_KEY uint64_t gva_key


static uint64_t const HASH_TABLE_NOT_FOUND = UINT64_MAX;


void*
hash_table_init(GVA_Allocator const allocator, size_t const capacity, size_t const item_size);


size_t
hash_table_index(void* const self, size_t const item_size, size_t const key);


void*
hash_table_ensure(GVA_Allocator const allocator, void* const self, size_t const item_size);


#define HASH_TABLE_DESTROY(allocator, self) ARRAY_DESTROY((allocator), (self))


#define HASH_TABLE_INDEX(self, key) hash_table_index((self), sizeof(*(self)), (key))


#define HASH_TABLE_SET(allocator, self, key, value) (                   \
    (self) = hash_table_ensure((allocator), (self), sizeof(*(self))),   \
    array_header(self)->length += 1,                                    \
    (self)[hash_table_index((self), sizeof(*(self)), (key))] = (value)  \
)


#endif  // GVA_HASH_TABLE_H
