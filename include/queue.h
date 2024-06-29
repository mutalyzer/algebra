#ifndef QUEUE_H
#define QUEUE_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>     // size_t
#include <stdint.h>     // uint32_t

#include "alloc.h"      // VA_Allocator


typedef struct VA_Queue VA_Queue;


VA_Queue*
va_queue_init(VA_Allocator const allocator, size_t const capacity);


VA_Queue*
va_queue_destroy(VA_Allocator const allocator, VA_Queue* const queue);


uint32_t
va_queue_dequeue(VA_Queue* const queue);


void
va_queue_enqueue(VA_Queue* const queue, uint32_t const value);


#ifdef __cplusplus
} // extern "C"
#endif

#endif
