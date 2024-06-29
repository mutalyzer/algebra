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
    VA_Queue* const queue = allocator.alloc(allocator.context, NULL, 0, sizeof(*queue) + capacity * sizeof(*queue->data));  // OVERFLOW
    if (queue == NULL)
    {
        return NULL;  // OOM
    } // if
    queue->head = 0;
    queue->tail = 0;
    queue->capacity = capacity;
    return queue;
} // va_queue_init


inline VA_Queue*
va_queue_destroy(VA_Allocator const allocator, VA_Queue* const queue)
{
    return allocator.alloc(allocator.context, queue, sizeof(*queue) + queue->capacity * sizeof(*queue->data), 0);
} // va_queue_destroy


inline uint32_t
va_queue_dequeue(VA_Queue* const queue)
{
    if (queue->head == queue->tail)
    {
        return -1;  // UNDERFLOW
    } // if
    uint32_t const value = queue->data[queue->head];
    queue->head = (queue->head + 1) % queue->capacity;
    return value;
} // va_queue_dequeue


inline void
va_queue_enqueue(VA_Queue* const queue, uint32_t const value)
{
    if ((queue->tail + 1) % queue->capacity == queue->head)
    {
        return;  // OVERFLOW
    } // if
    queue->data[queue->tail] = value;
    queue->tail = (queue->tail + 1) % queue->capacity;
} // va_queue_enqueue
