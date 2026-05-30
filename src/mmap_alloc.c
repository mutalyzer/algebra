// NOT FREESTANDING
// LINUX ONLY!
#define _GNU_SOURCE


#include <errno.h>      // DEBUG: errno
#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // DEBUG: perror

#include <fcntl.h>      // O_*, S_*, open
#include <sys/mman.h>   // MADV_, MAP_*, MREMAP_MAYMOVE, PROT_*,
                        // madvise, mmap, mremap, munmap
#include <sys/stat.h>   // stat, fstat
#include <unistd.h>     // ftruncate, close

#include "../include/mmap_alloc.h"      // GVA_Mmap_Context, gva_mmap_*
#include "common.h"     // MAX


inline void*
gva_mmap_allocate(void* const restrict context, void* const restrict ptr,
    size_t const old_size, size_t const new_size)
{
    (void) ptr;  // ignored (== ctx->addr)

    GVA_Mmap_Context* const restrict ctx = context;

    // free: NOP
    if (new_size == 0)
    {
        return NULL;
    } // if

    // deallocate or sufficient available: NOP
    if (new_size < old_size || new_size < ctx->len)
    {
        return ctx->addr;
    } // if

    // allocate more

    errno = 0;  // DEBUG
    if (ftruncate(ctx->fd, new_size) == -1)
    {
        perror("ftruncate()");  // DEBUG
        return NULL;  // error
    } // if

    errno = 0;  // DEBUG
    void* const restrict addr = mremap(ctx->addr, ctx->len, new_size, MREMAP_MAYMOVE);
    if (addr == MAP_FAILED)
    {
        perror("mremap()");  // DEBUG
        return NULL;  // error
    } // if
    ctx->addr = addr;
    ctx->len = new_size;

    return ctx->addr;
} // gva_mmap_allocate


GVA_Mmap_Context
gva_mmap_context_init(char const path[static 1], size_t const initial_len)
{
    static size_t const MINIMUM_LEN = 4096;

    GVA_Mmap_Context ctx = {NULL};

    errno = 0;  // DEBUG
    ctx.fd = open(path, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (ctx.fd == -1)
    {
        perror("open()");  // DEBUG
        return (GVA_Mmap_Context) {NULL};  // error
    } // if

    struct stat sb = {0};
    errno = 0;  // DEBUG
    if (fstat(ctx.fd, &sb) == -1)
    {
        perror("fstat()");  // DEBUG
        close(ctx.fd);  // unchecked
        return (GVA_Mmap_Context) {NULL};  // error
    } // if
    size_t const file_len = sb.st_size;

    ctx.len = MAX(MAX(MINIMUM_LEN, initial_len), file_len);

    if (file_len < ctx.len)
    {
        errno = 0;  // DEBUG
        if (ftruncate(ctx.fd, ctx.len) == -1)
        {
            perror("ftruncate()");  // DEBUG
            close(ctx.fd);  // unchecked
            return (GVA_Mmap_Context) {NULL};  // error
        } // if
    } // if

    errno = 0;  // DEBUG
    ctx.addr = mmap(NULL, ctx.len, PROT_READ | PROT_WRITE, MAP_POPULATE | MAP_SHARED, ctx.fd, 0);
    if (ctx.addr == MAP_FAILED)
    {
        perror("mmap()");  // DEBUG
        close(ctx.fd);  // unchecked
        return (GVA_Mmap_Context) {NULL};  // error
    } // if

    errno = 0;  // DEBUG
    if (madvise(ctx.addr, ctx.len, MADV_SEQUENTIAL | MADV_WILLNEED) == -1)
    {
        perror("madvise()");  // DEBUG
    } // if

    return ctx;
} // gva_mmap_context_init


inline void
gva_mmap_context_destroy(GVA_Mmap_Context const ctx[static 1])
{
    errno = 0;  // DEBUG
    if (munmap(ctx->addr, ctx->len) == -1)
    {
        perror("munmap()");  // DEBUG
    } // if
    errno = 0;  // DEBUG
    if (close(ctx->fd) == -1)
    {
        perror("close()");  // DEBUG
    } // if
} // gva_mmap_context_destroy
