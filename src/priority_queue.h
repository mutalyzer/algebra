#ifndef GVA_PRIORITY_QUEUE_H
#define GVA_PRIORITY_QUEUE_H


#include <stdbool.h>    // bool
#include <stddef.h>     // size_t

#include "../include/allocator.h"   // GVA_Allocator
#include "../include/types.h"       // gva_uint


typedef struct
{
    gva_uint idx;
    gva_uint included;
    gva_uint excluded;
} State;


typedef struct
{
    gva_uint*     heap;
    State*        states;
    size_t        capacity;
    GVA_Allocator allocator;
} Priority_Queue;


Priority_Queue
priority_queue_init(GVA_Allocator const allocator, size_t const capacity);


void
priority_queue_destroy(Priority_Queue self[static 1]);


bool
priority_queue_empty(Priority_Queue const self);


size_t
priority_queue_pop(Priority_Queue self[static 1]);


void
priority_queue_push(Priority_Queue self[static 1], size_t const key,
    size_t const included, size_t const excluded);


#endif  // GVA_PRIORITY_QUEUE_H
