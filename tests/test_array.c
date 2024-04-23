#include <assert.h>     // assert
#include <stddef.h>     // NULL, size_t
#include <stdint.h>     // uint32_t
#include <stdio.h>      // stderr, fprintf
#include <stdlib.h>     // EXIT_SUCCESS, free, realloc
#include <string.h>     // memcmp

#include "../include/alloc.h"       // VA_Allocator
#include "../include/std_alloc.h"   // std_allocator


#include "../src/array.c"


static inline void*
null_alloc(void* const context, void* const ptr, size_t const old_size, size_t const size)
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
        VA_Allocator const* const allocator;
        size_t const capacity;
        VA_Array const header;
    } const tests[] =
    {
        {&std_allocator,   0, {  4, 4, 0}},
        {&std_allocator,   1, {  4, 4, 0}},
        {&std_allocator,   4, {  4, 4, 0}},
        {&std_allocator,   5, {  5, 4, 0}},
        {&std_allocator, 100, {100, 4, 0}},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        uint32_t* array = va_array_init(tests[i].allocator, tests[i].capacity, sizeof(*array));
        assert(memcmp(((VA_Array*) array) - 1, &tests[i].header, sizeof(tests[i].header)) == 0);
        array = va_array_destroy(tests[i].allocator, array);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_array_init


static void
test_va_array_init_null(void)
{
    static struct {
        size_t const capacity;
    } const tests[] =
    {
        {1},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        assert(va_array_init(&null_allocator, tests[i].capacity, sizeof(uint32_t)) == NULL);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_array_init_null


static void
test_va_array_ensure(void)
{
    static struct {
        VA_Allocator const* const allocator;
        size_t const capacity;
        size_t const count;
        VA_Array const header;
    } const tests[] =
    {
        {&std_allocator,   0,   0, {  4, 4, 0}},
        {&std_allocator,   0,   1, {  4, 4, 0}},
        {&std_allocator,   0,   4, {  4, 4, 0}},
        {&std_allocator,   0,   5, {  8, 4, 0}},
        {&std_allocator, 100, 100, {100, 4, 0}},
        {&std_allocator, 100, 101, {200, 4, 0}},
        {&std_allocator,   0, 100, {128, 4, 0}},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        uint32_t* array = va_array_init(tests[i].allocator, tests[i].capacity, sizeof(*array));
        array = va_array_ensure(tests[i].allocator, array, tests[i].count);
        assert(memcmp(((VA_Array*) array) - 1, &tests[i].header, sizeof(tests[i].header)) == 0);
        array = va_array_destroy(tests[i].allocator, array);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_array_ensure


static void
test_va_array_ensure_null(void)
{
    static struct {
        size_t const capacity;
        size_t const count;
    } const tests[] =
    {
        {0, 5},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        uint32_t* array = va_array_init(&std_allocator, tests[i].capacity, sizeof(*array));
        assert(va_array_ensure(&null_allocator, array, tests[i].count) == NULL);
        array = va_array_destroy(&std_allocator, array);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_array_ensure_null


static void
test_va_array_destroy(void)
{
    static struct {
        VA_Allocator const* const allocator;
        size_t const capacity;
    } const tests[] =
    {
        {&std_allocator, 0},
        {&null_allocator, 0},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        uint32_t* array = va_array_init(tests[i].allocator, tests[i].capacity, sizeof(*array));
        array = va_array_destroy(&std_allocator, array);
        assert(array == NULL);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_array_destroy


static void
test_va_array_append(void)
{
    static struct {
        VA_Allocator const* const allocator;
        size_t const capacity;
        uint32_t const count;
        VA_Array const header;
    } const tests[] =
    {
        {&std_allocator, 0,  0, { 4, 4,  0}},
        {&std_allocator, 0, 10, {16, 4, 10}},
        {&std_allocator, 5, 10, {10, 4, 10}},
    }; // tests

    fprintf(stderr, "%s:%s: ", __FILE__, __func__);
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        uint32_t* array = va_array_init(tests[i].allocator, tests[i].capacity, sizeof(*array));
        for (uint32_t j = 0; j < tests[i].count; ++j)
        {
            va_array_append(tests[i].allocator, array, j);
        } // for
        assert(memcmp(((VA_Array*) array) - 1, &tests[i].header, sizeof(tests[i].header)) == 0);
        for (uint32_t j = 0; j < tests[i].count; ++j)
        {
            assert(array[j] == j);
        } // for
        array = va_array_destroy(&std_allocator, array);
        fprintf(stderr, ".");
    } // for
    fprintf(stderr, "  passed\n");
} // test_va_array_append


int
main(void)
{
    test_va_array_init();
    test_va_array_init_null();
    test_va_array_ensure();
    test_va_array_ensure_null();
    test_va_array_destroy();
    test_va_array_append();
    return EXIT_SUCCESS;
} // main
