#include <stdbool.h>    // bool
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint32_t

#include "../include/alloc.h"       // VA_Allocator
#include "../include/queue.h"       // va_queue_*


struct VA_Queue
{
    size_t head;
    size_t tail;
    size_t capacity;
    uint32_t data[];
};


inline VA_Queue*
va_queue_init(VA_Allocator const allocator, size_t const capacity)
{
    VA_Queue* const self = allocator.alloc(allocator.context, NULL, 0, sizeof(*self) + capacity * sizeof(*self->data));  // OVERFLOW
    if (self == NULL)
    {
        return NULL;  // OOM
    } // if
    self->head = 0;
    self->tail = 0;
    self->capacity = capacity;
    return self;
} // va_queue_init


inline VA_Queue*
va_queue_destroy(VA_Allocator const allocator, VA_Queue* const self)
{
    return allocator.alloc(allocator.context, self, sizeof(*self) + self->capacity * sizeof(self->data[0]), 0);
} // va_queue_destroy


inline uint32_t
va_queue_dequeue(VA_Queue* const self)
{
    if (self->head == self->tail)
    {
        return -1;  // UNDERFLOW
    } // if
    uint32_t const value = self->data[self->head];
    self->head = (self->head + 1) % self->capacity;
    return value;
} // va_queue_dequeue


inline bool
va_queue_enqueue(VA_Queue* const self, uint32_t const value)
{
    if ((self->tail + 1) % self->capacity == self->head)
    {
        return false;  // OVERFLOW
    } // if
    self->data[self->tail] = value;
    self->tail = (self->tail + 1) % self->capacity;
    return true;
} // va_queue_enqueue


inline bool
va_queue_is_empty(VA_Queue const* const self)
{
    return self->head == self->tail;
} // va_queue_is_empty
