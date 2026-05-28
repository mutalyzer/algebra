// NOT FREESTANDING
#ifndef GVA_MMAP_H
#define GVA_MMAP_H


// FIXME: proper header + source


#include <errno.h>      // errno
#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // stderr, stdout, fprintf
#include <stdlib.h>     // EXIT_*, free, malloc
#include <string.h>     // strerror

#include <fcntl.h>      // O_*, S_*, open
#include <sys/mman.h>   // MADV_, MAP_*, MREMAP_MAYMOVE, MS_ASYNC, PROT_*,
                        // madvise, mmap, mremap, msync
#include <sys/stat.h>   // stat, fstat
#include <unistd.h>     // ftruncate, close

#include "../include/allocator.h"      // GVA_Allocator
#include "array.h"      // ARRAY*, array_length


typedef struct
{
    void*  addr;
    size_t len;
    int    fd;
} MMAP_Context;


static inline void*
mmap_allocate(void* const restrict context, void* const restrict ptr, size_t const old_size, size_t const new_size)
{
    (void) ptr;  // == ctx->addr

    MMAP_Context* const ctx = context;
    ctx->len = new_size;

    if (new_size == 0)
    {
        return NULL;
    } // if

    if (ctx->addr == NULL)
    {
        errno = 0;
        ctx->addr = mmap(NULL, ctx->len, PROT_READ | PROT_WRITE, MAP_POPULATE | MAP_SHARED, ctx->fd, 0);
        if (ctx->addr == MAP_FAILED)
        {
            fprintf(stderr, "mmap(): %s\n", strerror(errno));
            return NULL;
        } // if
    } // if
    else
    {
        errno = 0;
        ctx->addr = mremap(ctx->addr, old_size, ctx->len, MREMAP_MAYMOVE);
        if (ctx->addr == MAP_FAILED)
        {
            fprintf(stderr, "mremap(): %s\n", strerror(errno));
            return NULL;
        } // if
    } // else

    errno = 0;
    if (ftruncate(ctx->fd, ctx->len))
    {
        fprintf(stderr, "ftruncate(): %s\n", strerror(errno));
        return NULL;
    } // if

    return ctx->addr;
} // mmap_allocate


static GVA_Allocator
mmap_allocator_init(char const path[static 1])
{
    MMAP_Context* const ctx = malloc(sizeof(*ctx));  // FIXME
    if (ctx == NULL)
    {
        return (GVA_Allocator) {NULL};  // OOM
    } // if
    ctx->addr = NULL;
    ctx->len = 0;

    errno = 0;
    ctx->fd = open(path, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (ctx->fd == -1)
    {
        fprintf(stderr, "open(): %s\n", strerror(errno));
        free(ctx);
        return (GVA_Allocator) {NULL};
    } // if

    struct stat sb = {0};
    errno = 0;
    if (fstat(ctx->fd, &sb) == -1)
    {
        fprintf(stderr, "fstat(): %s\n", strerror(errno));
        close(ctx->fd);
        free(ctx);
        return (GVA_Allocator) {NULL};
    } // if

    if (sb.st_size > 0)
    {
        ctx->len = sb.st_size;

        errno = 0;
        ctx->addr = mmap(NULL, ctx->len, PROT_READ | PROT_WRITE, MAP_POPULATE | MAP_SHARED, ctx->fd, 0);
        if (ctx->addr == MAP_FAILED)
        {
            fprintf(stderr, "mmap(): %s\n", strerror(errno));
            close(ctx->fd);
            free(ctx);
            return (GVA_Allocator) {NULL};
        } // if

        errno = 0;
        if (madvise(ctx->addr, ctx->len, MADV_SEQUENTIAL | MADV_WILLNEED))
        {
            // not critical
            fprintf(stderr, "madvise: %s\n", strerror(errno));
        } // if
    } // if

    return (GVA_Allocator)
    {
        .allocate = mmap_allocate,
        .context = ctx,
    };
} // mmap_allocator_init


static inline void
mmap_allocator_destroy(GVA_Allocator const allocator)
{
    MMAP_Context* const ctx = allocator.context;
    errno = 0;
    if (msync(ctx->addr, ctx->len, MS_ASYNC))
    {
        fprintf(stderr, "msync(): %s\n", strerror(errno));
    } // if
    close(ctx->fd);
    free(allocator.context);  // FIXME
} // mmap_allocator_destroy


static inline void*
array_load(void* const ptr)
{
    if (ptr == NULL)
    {
        return NULL;
    } // if

    MMAP_Context* const ctx = ptr;
    return ctx->addr;
} // array_load

#endif  // GVA_MMAP_H
