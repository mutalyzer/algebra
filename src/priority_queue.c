#include <stdbool.h>    // bool
#include <stddef.h>     // size_t

#include "../include/allocator.h"   //  GVA_Allocator
#include "../include/types.h"       // GVA_NULL, gva_uint

#include "array.h"              // ARRAY_*
#include "priority_queue.h"     // Priority_Queue, priority_queue_*, State


static inline bool
greater_than(size_t const lhs_included, size_t const lhs_excluded,
    size_t const rhs_included, size_t const rhs_excluded)
{
    return lhs_included > rhs_included ||
        (lhs_included == rhs_included && lhs_excluded <= rhs_excluded);
} // greater_than


static inline void
swap(gva_uint lhs[static restrict 1], gva_uint rhs[static restrict 1])
{
    gva_uint const temp = *lhs;
    *lhs = *rhs;
    *rhs = temp;
} // swap


inline Priority_Queue
priority_queue_init(GVA_Allocator const allocator, size_t const capacity)
{
    State* const states = allocator.allocate(allocator.context, NULL, 0, capacity * sizeof(*states));
    if (states == NULL)
    {
        return (Priority_Queue) {NULL};
    } // if

    for (size_t i = 0; i < capacity; ++i)
    {
        states[i].idx = GVA_NULL;
    } // for

    return (Priority_Queue)
    {
        .states = states,
        .capacity = capacity,
        .allocator = allocator,
    };
} // priority_queue_init


inline void
priority_queue_destroy(Priority_Queue self[static 1])
{
    self->heap = ARRAY_DESTROY(self->allocator, self->heap);
    self->states = self->allocator.allocate(self->allocator.context, self->states, self->capacity * sizeof(State), 0);
    self->capacity = 0;
} // priority_queue_destroy


inline bool
priority_queue_empty(Priority_Queue const self)
{
    return array_length(self.heap) == 0;
} // priority_queue_empty


inline size_t
priority_queue_pop(Priority_Queue self[static 1])
{
    size_t const idx = self->heap[0];
    self->states[self->heap[0]].idx = GVA_NULL;
    size_t const len = array_header(self->heap)->length -= 1;
    self->heap[0] = self->heap[len];
    if (len > 0)
    {
        self->states[self->heap[0]].idx = 0;
    } // if

    size_t i = 0;
    while (i < len / 2)
    {
        size_t child = 2 * i + 1;  // left child
        if (child + 1 < len && greater_than(self->states[self->heap[child + 1]].included, self->states[self->heap[child + 1]].excluded,
                                            self->states[self->heap[child]].included, self->states[self->heap[child]].excluded))
        {
            child += 1;  // right child
        } // if
        if (greater_than(self->states[self->heap[i]].included, self->states[self->heap[i]].excluded,
                         self->states[self->heap[child]].included, self->states[self->heap[child]].excluded))
        {
            break;
        } // if
        swap(&self->states[self->heap[i]].idx, &self->states[self->heap[child]].idx);
        swap(&self->heap[i], &self->heap[child]);
        i = child;
    } // while

    return idx;
} // priority_queue_pop


inline void
priority_queue_push(Priority_Queue self[static 1], size_t const key,
    size_t const included, size_t const excluded)
{
    if (self->states[key].idx == GVA_NULL)
    {
        self->states[key] = (State)
        {
            .idx = ARRAY_APPEND(self->allocator, self->heap, key),
            .included = included,
            .excluded = excluded,
        };
    } // if
    else
    {
        if (greater_than(self->states[key].included, self->states[key].excluded, included, excluded))
        {
            return;
        } // if
        self->states[key].included = included;
        self->states[key].excluded = excluded;
    } // else

    size_t i = self->states[key].idx;
    while (i > 0)
    {
        size_t const parent = (i - 1) / 2;
        if (greater_than(self->states[self->heap[parent]].included, self->states[self->heap[parent]].excluded, included, excluded))
        {
            break;
        } // if
        swap(&self->states[self->heap[i]].idx, &self->states[self->heap[parent]].idx);
        swap(&self->heap[i], &self->heap[parent]);
        i = parent;
    } // while
} // priority_queue_push
