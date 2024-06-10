#include <assert.h>     // assert
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint32_t
#include <stdio.h>      // stderr, fprintf
#include <stdlib.h>     // EXIT_SUCCESS
#include <string.h>     // memcmp

#include "../include/alloc.h"       // VA_Allocator
#include "../include/std_alloc.h"   // va_std_allocator


#include "../src/array.c"


static inline void*
null_alloc(void* const restrict context, void* const restrict ptr, size_t const old_size, size_t const size)
{
    (void) context;
    (void) ptr;
    (void) old_size;
    (void) size;
    return NULL;
} // null_alloc


static VA_Allocator const null_allocator = { .alloc = null_alloc };


static void
test_va_array_init(void)
{
    static struct {
        size_t const capacity;
        VA_Array const header;
    } const tests[] =
    {
        {  1, {  1, 0}},
        { 10, { 10, 0}},
        {100, {100, 0}},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        uint32_t* array = va_array_init(va_std_allocator, tests[i].capacity, sizeof(*array));
        assert(memcmp(va_array_header(array), &tests[i].header, sizeof(tests[i].header)) == 0);
        array = va_array_destroy(va_std_allocator, array);
        assert(array == NULL);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_array_init


static void
test_va_array_init_zero(void)
{
    static struct {
        VA_Allocator const allocator;
        size_t const capacity;
    } const tests[] =
    {
        {va_std_allocator, 0},
        {null_allocator,   0},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        uint32_t* array = va_array_init(tests[i].allocator, tests[i].capacity, sizeof(*array));
        assert(array == NULL);
        array = va_array_destroy(tests[i].allocator, array);
        assert(array == NULL);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_array_init_zero


static void
test_va_array_init_null(void)
{
    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    uint32_t* array = va_array_init(null_allocator, 100, sizeof(*array));
    assert(array == NULL);
    array = va_array_destroy(null_allocator, array);
    assert(array == NULL);
    fprintf(stderr, ".");
    fprintf(stderr, "  passed\n");
} // test_va_array_init_null


static void
test_va_array_append(void)
{
    uint32_t const n = 100;
    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    uint32_t* array = NULL;
    for (uint32_t i = 0; i < n; ++i)
    {
        va_array_append(va_std_allocator, array, i);
    } // for
    assert(va_array_length(array) == n);
    for (uint32_t i = 0; i < va_array_length(array); ++i)
    {
        assert(array[i] == i);
    } // for
    array = va_array_destroy(va_std_allocator, array);
    assert(array == NULL);
    fprintf(stderr, ".");
    fprintf(stderr, "  passed\n");
} // test_va_array_append


int
main(void)
{
    test_va_array_init();
    test_va_array_init_zero();
    test_va_array_init_null();
    test_va_array_append();
    return EXIT_SUCCESS;
} // main
