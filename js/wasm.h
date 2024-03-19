#ifndef WASM_H
#define WASM_H


#include <stddef.h>     // size_t


#define WASM_EXPORT(name) __attribute__((export_name(name)))


extern void  console_log(size_t, char const*);

extern void  free(void*);
extern void* malloc(size_t);


WASM_EXPORT("__word_size")
size_t
__word_size(void)
{
    return sizeof(size_t);
} // __word_size


void*
calloc(size_t count, size_t size)
{
    return malloc(count * size);
} // calloc


__attribute__((noinline))
void*
realloc(void* ptr, size_t size)
{
    free(ptr);
    return malloc(size);
} // realloc


void*
memcpy(void* dest, void const* src, size_t n)
{
    char* bdest = dest;
    char const* bsrc = src;

    for (size_t i = 0; i < n; ++i)
    {
        bdest[i] = bsrc[i];
    } // for
    return dest;
} // memcpy


#endif
